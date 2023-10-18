"""
Scripts for analyzing the results of Wright-Fisher evolutionary simulations. 
"""

import numpy as np
import pandas as pd
from tqdm.auto import tqdm


def get_genotype_frequencies(generations,cutoff=50):
    """
    Extract genotype frequencies over simulation. 
    
    Parameters
    ----------
    generations : list
        list of dicts holding results returned by a Wright-Fisher simulation
    cutoff : int, default=50
        drop genotypes seen fewer than cutoff times over the simulation
    
    Returns
    -------
    df : pandas.DataFrame
        dataframe with genotype frequencies over the course of the simulation. 
        columns are genotypes; rows are generations. 
    """

    max_counts = {}
    population_size = None
    for g in tqdm(generations):

        if population_size is None:
            population_size = 0
            for k in g:
                population_size += g[k]
                
        for k in g:
            if k not in list(max_counts.keys()):
                max_counts[k] = g[k]
                
            if g[k] > max_counts[k]:
                max_counts[k] = g[k]
                
    keys_to_keep = []
    for k in max_counts:
        if max_counts[k] > cutoff:
            keys_to_keep.append(k)
    
    out_dict = {}
    for k in keys_to_keep:
        out_dict[k] = np.zeros(len(generations),dtype=float)
        for i in range(len(generations)):
            try:
                out_dict[k][i] = generations[i][k]/population_size
            except KeyError:
                pass
        
    return pd.DataFrame(out_dict)
    

def get_most_common(some_dict):
    """
    Return key/value pair where value is highest. 
    """
    
    keys = list(some_dict.keys())
    values = [some_dict[k] for k in keys]
    to_sort = list(zip(values,keys))
    to_sort.sort()
    
    return to_sort[-1][1], to_sort[-1][0]

def extract_alignment(gc,tree,include_ancestors=False):
    """
    Extract an alignment from a simulation done on an evolutionary tree.
    
    Parameters
    ----------
    gc : GenotypeContainer
        genotype container with all genotypes seen in the simulation
    tree : ete3.Tree
        tree object with node.generations features (output by the simulation)
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

    # Get wildtype sequence
    wt_seq = list(gc.wt_sequence)

    out = {}
    counter = 0
    
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
            for m in gc.genotypes[genotype].mutations:
                idx = int(m[1:-1]) - 1
                aa = m[-1]
                this_seq[idx] = aa

            # Update output
            out[n.name] = "".join(this_seq)

    return out, tree