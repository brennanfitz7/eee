
import eee
from eee.calcs import read_json
from eee.analysis.wf.analysis import get_most_common
from eee._private.check.standard import check_bool

import os
import glob
import pickle

def load_wf_tree_sim(calc_dir,
                     keep_all_generations=False):
    """
    Load the output from a WrightFisherTreeSimulation calculation.
    
    Parameters
    ----------
    calc_dir : str
        directory holding calculation
    keep_all_generations : bool, default=False
        keep every generation along every branch. If False, store only the 
        last generation at each node. (Warning: setting to True can lead to a 
        huge tree in memory). 
        
    Returns
    -------
    sc : SimulationContainer subclass
        initialized SimulationContainer subclass with ensemble, fitness, and
        ddg loaded.
    sc_kwargs : dict
        dictionary with run parameters. sc.run(**calc_params) will run the 
        calculation.
    tree : ete3.Tree
        tree with calculation generations loaded into the "generations" 
        node feature. 
    genotypes : pandas.DataFrame
        dataframe with genotypes seen over the simulations. 
    """

    # Check argument stanity
    calc_dir = f"{calc_dir}"
    if not os.path.isdir(calc_dir):
        err = f"\ncalc_dir ({calc_dir}) should be a simulation output directory\n\n"
        raise ValueError(err)
    
    keep_all_generations = check_bool(keep_all_generations,
                                      variable_name="keep_all_generations")


    # Load calculation
    sc, sc_kwargs = read_json(os.path.join(calc_dir,
                                           "input",
                                           "simulation.json"),
                              tree_fmt=3)
    base = sc_kwargs["write_prefix"] 
    tree = sc_kwargs["tree"]

    # Figure out name of the root ancestor (kind of a hack because ete3 does 
    # not write root name to newick). Pull from name of burn-in pickle file.
    burn_in = glob.glob(os.path.join(calc_dir,
                                    f"{base}*burn-in*.pickle"))[0]
    anc_root = burn_in.split("burn-in-")[1].split(".pickle")[0]

    # Go through tree
    for n in tree.traverse():
        
        # If root, give it the appropriate name
        if n.is_root():
            n.name = anc_root
        
        # Get ancestors
        ancestors = n.get_ancestors()
        
        # No ancestors - this is root; get burn-in pickle
        if len(ancestors) == 0:
            branch = f"burn-in-{n.name}"
        
        # Get ancestor to get appropriate pickle file name
        else:
            branch = f"{ancestors[0].name}-{n.name}"
        
        # Load generations from pickle
        pickle_file = os.path.join(calc_dir,f"{base}_{branch}.pickle")
        with open(pickle_file,"rb") as f:
            generations = pickle.load(f)

        # Drop all generations but the last one (corresponding to this node) 
        if not keep_all_generations:
            generations = [generations[-1]]

        # append generations to pickle
        n.add_feature("generations",generations)

    # Load the genotypes spreadsheet
    genotypes = eee.io.read_dataframe(os.path.join(calc_dir,f"{base}_genotypes.csv"))

    return sc, sc_kwargs, tree, genotypes


def extract_alignment(calc_dir,
                      include_ancestors=False):
    """
    Extract an alignment from a simulation done on an evolutionary tree.
    
    Parameters
    ----------
    calc_dir : str
        directory holding the calculation
    include_ancestors : bool, default=False
        whether or not to include ancestors in the alignment output. If True,
        ancestors are given numbered names of the form anc{counter}. The 
        ancestral nodes in the tree are given node.name values with the 
        appropriate name. 
    
    Returns
    -------
    out : dict
        dictionary where keys are node names and values are sequences (as str)
    tree : ete3.Tree
        tree object used for the inference. If include_ancestors = True, this 
        will have added ancestor names. 
    """

    sc, sc_kwargs, tree , genotypes = load_wf_tree_sim(calc_dir)
    
    # Get wildtype sequence
    wt_seq = list(sc.gc.wt_sequence)

    out = {}
    counter = 0
    
    genotype_dict = dict(zip(genotypes["genotype"],
                             genotypes["mutations"]))

    # Go through tree
    for n in tree.traverse():

        # If a leaf (or we are including ancestors)...
        if n.is_leaf() or include_ancestors:

            # Give ancestor name to unlabeled nodes
            if n.name == "":
                n.name = f"anc{counter}"
                counter += 1

            # Start with wildtype sequence
            this_seq = wt_seq[:]

            # Get most common genotype in the last generation
            genotype, _ = get_most_common(n.generations[-1])
            
            # Update this_seq with the mutations present in that genotype
            muts = genotype_dict[genotype]
            if issubclass(type(muts),str):
                muts = muts.split("/")
            else:
                # wt -- is nan
                muts = []

            for m in muts:
                idx = int(m[1:-1]) - 1
                aa = m[-1]
                this_seq[idx] = aa

            # Update output
            out[n.name] = "".join(this_seq)

    out_fasta = []
    for k in out:
        out_fasta.append(f">{k}\n{out[k]}")

    out_fasta = "\n".join(out_fasta)

    return out, out_fasta