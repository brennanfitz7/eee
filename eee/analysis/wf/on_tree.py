
from eee.analysis.wf.analysis import get_most_common

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
    pass

    # # Get wildtype sequence
    # wt_seq = list(gc.wt_sequence)

    # out = {}
    
    # # Go through tree
    # for n in tree.traverse():

    #     # If a leaf (or we are including ancestors)...
    #     if n.is_leaf() or include_ancestors:

    #         # Start with wildtype sequence
    #         this_seq = wt_seq[:]



    #         # Get most common genotype in the last generation
    #         genotype, _ = get_most_common(n.population)

    #         # Update this_seq with the mutations present in that genotype
    #         for m in gc.genotypes[genotype].mutations:
    #             idx = int(m[1:-1]) - 1
    #             aa = m[-1]
    #             this_seq[idx] = aa

    #         # Update output
    #         out[n.name] = "".join(this_seq)

    # return out, tree