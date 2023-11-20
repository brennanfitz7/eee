"""
Simulate evolution of an ensemble along an evolutionary tree. 
"""

from eee.core.genotype import Genotype
from eee.core.engine import wright_fisher
from eee.core.engine.wright_fisher import get_num_accumulated_mutations
from eee.core.engine.wright_fisher import write_wf_outputs

from eee._private.check.eee import check_num_generations
from eee._private.check.eee import check_mutation_rate
from eee._private.check.eee import check_burn_in_generations
from eee._private.check.eee import check_wf_population

import numpy as np
from tqdm.auto import tqdm
import ete3

import pickle

def _simulate_branch(start_node,
                     end_node,
                     gc,
                     mutation_rate,
                     num_generations,
                     write_prefix,
                     rng):
    """
    Simulate evolution along a branch. start_node and end_node are ete3.Tree 
    nodes. Simulate from one to the other, storing the final generation as 
    the population feature in the end-node. Write out a pickle file with 
    the entire simulation along the branch. 
    """
    
    starting_pop = start_node.population

    # Get number of mutations to accumulate based on the branch length times 
    # the sequence length
    branch_length = start_node.get_distance(end_node)
    sequence_length = len(gc.wt_sequence)
    num_mutations = int(np.round(branch_length*sequence_length,0))

    # Get the number of mutations at the starting node. The number to accumulate
    # over the branch is the start + the branch length. (This is the total number
    # of mutations that have accumulated, including reversions and multiple 
    # mutations at the same site). 
    num_mutations_start = get_num_accumulated_mutations(seen=list(starting_pop.keys()),
                                                        counts=list(starting_pop.values()),
                                                        gc=gc)
    
    num_mutations = num_mutations + num_mutations_start

    # Force at least one mutation to occur on the branch
    if num_mutations == 0:
        num_mutations = 1
    
    gc, generations = wright_fisher(gc,
                                    population=starting_pop,
                                    mutation_rate=mutation_rate,
                                    num_generations=num_generations,
                                    num_mutations=num_mutations,
                                    verbose=False,
                                    write_prefix=None,
                                    rng=rng)
    
    # record generations. (Note, in pickle files, ) 
    pickle_name = f"{write_prefix}_{start_node.name}-{end_node.name}.pickle"
    with open(pickle_name,"wb") as f:
        pickle.dump(generations,f)

    end_node.add_feature("population",generations[-1])


def follow_tree(gc,
                tree,
                population=1000,
                mutation_rate=0.01,
                num_generations=100,
                burn_in_generations=10,
                write_prefix="eee_follow-tree",
                rng=None):
    """
    Run a Wright-Fisher simulation following an evolutionary tree. 
    
    Parameters
    ----------
    gc : GenotypeContainer
        container holding genotypes with either a wildtype sequence or a 
        sequences from a previous evolutionary simulation. This will bring in 
        everything we need to calculate fitness of each genotype. 
    tree : ete3.Tree
        ete3 tree with branch lengths and tip labels
    population : dict, list-like, or int
        population for the simulation. Can be a dictionary of populations where
        keys are genotypes and values are populations. Can be an array where the
        length is the total population size and the values are the genotypes. 
        {5:2,4:1,0:2} is equivalent to [5,5,4,0,0]. If an integer, interpret as
        the population size. Create an initial population consisting of all 
        wildtype.
    mutation_rate : float
        Lambda for poisson distribution to select the number of genotypes to 
        mutate at each generation. Should be >= 0. 
    num_generations : int
        number of generations to run the simulation. Should be >= 1. If 
        num_mutations is specified, this is the *maximum* number of generations 
        allowed. 
    num_mutations : int, optional
        stop the simulation after the most frequent genotype has num_mutations 
        mutations. Should be >= 1. If specified, the simulation will run until
        either num_mutations is reached OR the simulation hits num_generations.
    burn_in_generations : int, default=10
        run a Wright-Fisher simulation burn_in_generations long to generate an
        ancestral population. Must be >= 0. 
    write_prefix : str
        write output files during the run with this prefix. 
    rng : numpy.random._generator.Generator, optional
        random number generator object to allow reproducible sims. If None, one
        is created locally. 

    Returns
    -------
    gc : GenotypeContainer
        updated GenotypeContainer holding all new genotypes that appeared over 
        the simulation
    generations : list
        list of dicts where each entry in the list is a generation. each dict
        holds the population at that generation with keys as genotypes and 
        values as populations. 
    
    """

    if not issubclass(type(gc),Genotype):
        err = "\ngc must be of type Genotype\n\n"
        raise ValueError(err)    
    
    if not issubclass(type(tree),ete3.TreeNode):
        err = "\ntree must be an ete3.Tree object\n\n"
        raise ValueError(err)  

    population = check_wf_population(population)
    population_size = len(population)
    mutation_rate = check_mutation_rate(mutation_rate)
    num_generations = check_num_generations(num_generations)
    burn_in_generations = check_burn_in_generations(burn_in_generations)
    
    if write_prefix is None:
        err = "\nwrite_prefix must not be None\n\n"
        raise ValueError(err)
    write_prefix = f"{write_prefix}"

    if rng is None:
        rng = np.random.Generator(np.random.PCG64())

    # Figure out the number of branches for the status bar
    total_branches = 1
    num_ancestors = 0
    for n in tree.traverse(strategy="levelorder"):
        if not n.is_leaf():
            total_branches += 2
            num_ancestors += 1

    # Get format string for ancestor names
    num_positions = len(f"{num_ancestors}") + 1
    anc_fmt_string = "anc{0:" + f"0{num_positions}" + "d}"

    gc_filename = f"{write_prefix}_genotypes.csv"

    genotypes_to_keep = []

    pbar = tqdm(total=total_branches)
    with pbar:

        # burn in to generate initial population
        gc, generations = wright_fisher(gc=gc,
                                        population=population_size,
                                        mutation_rate=mutation_rate,
                                        num_generations=burn_in_generations,
                                        verbose=False,
                                        write_prefix=None,
                                        rng=rng)

        pbar.update(n=1)

        # record generations. (Note, the population of the ancestral node and
        # the derived node is stored in each pickle file -- the same generation
        # is thus repeated across pickle files). 
        anc_name = anc_fmt_string.format(0)
        pickle_name = f"{write_prefix}_burn-in-{anc_name}.pickle"
        with open(pickle_name,"wb") as f:
            pickle.dump(generations,f)

        anc_counter = 0

        # Get the tree root and append the generations from the burn in.
        root = tree.get_tree_root()
        root.add_feature("population",generations[-1])
        root.name = anc_fmt_string.format(anc_counter)
        anc_counter += 1

        # Dump genotypes
        genotypes_to_keep = list(generations[-1].keys())
        gc.dump_to_csv(filename=gc_filename,
                       keep_genotypes=genotypes_to_keep)

        for n in tree.traverse(strategy="levelorder"):
            
            if not n.is_leaf():

                # Get descendants
                left, right = n.get_children()

                # Simulate evolution from n to left descendant. 
                if left.name == "":
                    left.name = anc_fmt_string.format(anc_counter)
                    anc_counter += 1

                _simulate_branch(start_node=n,
                                 end_node=left,
                                 gc=gc,
                                 mutation_rate=mutation_rate,
                                 num_generations=num_generations,
                                 write_prefix=write_prefix,
                                 rng=rng)
                
                # Dump genotypes to file
                genotypes_to_keep.extend(list(left.population.keys()))
                genotypes_to_keep = list(set(genotypes_to_keep))
                gc.dump_to_csv(filename=gc_filename,
                               keep_genotypes=genotypes_to_keep)

                pbar.update(n=1)

                # Simulate evolution from n to right descendent. (implicitly updates
                # gc and right node)
                if right.name == "":
                    right.name = anc_fmt_string.format(anc_counter)
                    anc_counter += 1

                _simulate_branch(start_node=n,
                                 end_node=right,
                                 gc=gc,
                                 mutation_rate=mutation_rate,
                                 num_generations=num_generations,
                                 write_prefix=write_prefix,
                                 rng=rng)
                
                # Dump genotypes to file
                genotypes_to_keep.extend(list(right.population.keys()))
                genotypes_to_keep = list(set(genotypes_to_keep))
                gc.dump_to_csv(filename=gc_filename,
                               keep_genotypes=genotypes_to_keep)

                pbar.update(n=1)
                

    # Write tree
    tree.write(format=3,
               outfile=f"{write_prefix}.newick")

    # Write genotypes
    gc.dump_to_csv(filename=gc_filename,
                   keep_genotypes=None)

    return gc, tree