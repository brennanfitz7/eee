
import eee

from eee.calcs import read_json
from eee.analysis.wf.analysis import get_most_common
from eee.analysis.wf.analysis import read_genotypes_file
from eee._private.check.standard import check_bool

from eee.core.data import AA_TO_INT
from eee.core.data import INT_TO_AA

import numpy as np
import pandas as pd
from tqdm.auto import tqdm

import os
import glob
import pickle
import itertools


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
        last generation at each node. (Warning: setting to True can lead to 
        crashes due to huge numbers of genotypes in memory). 
        
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

    # Check argument sanity
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

    genotypes_seen = set()

    total_branches = 1
    for n in tree.traverse(strategy="levelorder"):
        if not n.is_leaf():
            total_branches += 2

    print(f"Reading branch trajectories",flush=True)

    pbar = tqdm(total=total_branches)
    with pbar:

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
                genotypes = set(generations[0].keys())
                genotypes_seen = genotypes_seen.union(genotypes)
            
            # append generations to pickle
            n.add_feature("generations",generations)

            pbar.update(n=1)

    # Load the genotypes spreadsheet
    genotypes_file = os.path.join(calc_dir,f"{base}_genotypes.csv")
    if keep_all_generations:
        genotypes_seen = None
    genotypes = read_genotypes_file(genotypes_file=genotypes_file,
                                    only_genotypes=genotypes_seen)

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

    sc, sc_kwargs, tree , genotypes = load_wf_tree_sim(calc_dir,
                                                       keep_all_generations=False)
    
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

def alignment_mic(ali_dict):
    """
    Calculate mutual information (log2) between columns in an alignment.
    
    Parameters
    ----------
    ali_dict : dict
        alignment dictionary where values are aligned sequences as strings
    
    Returns
    -------
    main_df : pandas.DataFrame
        dataframe with entry for each i,j pair of columns/amino acids in the 
        alignment
    mic_df : pandas.DataFrame
        dataframe with the mic for each i,j pair of columns in the alignment
    """

    num_seqs = len(ali_dict)
    num_cols = len(ali_dict[list(ali_dict.keys())[0]])

    # Load alignment into a numpy array
    alignment = np.zeros((num_seqs,num_cols),dtype=int)
    for i, a in enumerate(ali_dict):
        alignment[i,:] = [AA_TO_INT[c] for c in list(ali_dict[a])]

    # Create dictionaries to store outputs (for ultimate conversion to a dataframe)
    cols = ["si","sj","aai","aaj","Ni","Nj","Nij","Pi","Pj","Pij"]
    main_out = dict([(c,[]) for c in cols])
    mic_out = {"si":[],"sj":[],"mic":[]}

    # Iterate over pairs of sites i,j
    for i in range(num_cols):

        # Create dictionary of amino acid counts for column i
        values, counts = np.unique(alignment[:,i],return_counts=True)
        i_singles_seen = dict(zip(values,counts))

        for j in range(i+1,num_cols):

            # Create dictionary of amino acid counts for column j
            values, counts = np.unique(alignment[:,j],return_counts=True)
            j_singles_seen = dict(zip(values,counts))
            
            # Create dictionary of paired amino acid counts at columns i,j
            pairs, counts = np.unique(alignment[:,[i,j]],axis=0,return_counts=True)
            ij_pairs_seen = dict(zip([tuple(p) for p in pairs],counts))

            # Calculate mic for each amino acid pair seen
            Pi_vec = np.zeros(len(ij_pairs_seen))
            Pj_vec = np.zeros(len(ij_pairs_seen))
            Pij_vec = np.zeros(len(ij_pairs_seen))
            for pair_counter, pair in enumerate(ij_pairs_seen):

                # Counts
                n_i = i_singles_seen[pair[0]]
                n_j = j_singles_seen[pair[1]]
                n_ij = ij_pairs_seen[pair]

                # probabilities --> P_ij is P(ij|j)*P(j)
                P_i = n_i/num_seqs
                P_j = n_j/num_seqs
                P_ij = n_ij/n_j*P_j
                
                # Record results
                main_out["si"].append(i)
                main_out["sj"].append(j)
                main_out["aai"].append(INT_TO_AA[pair[0]])
                main_out["aaj"].append(INT_TO_AA[pair[1]])
                main_out["Ni"].append(n_i)
                main_out["Nj"].append(n_j)
                main_out["Nij"].append(n_ij)
                main_out["Pi"].append(P_i)
                main_out["Pj"].append(P_j)
                main_out["Pij"].append(P_ij)

                # Vector of Ps for mic calculation
                Pi_vec[pair_counter] = P_i
                Pj_vec[pair_counter] = P_j
                Pij_vec[pair_counter] = P_ij

            # Append mic for whole site. 
            mic = np.sum(Pij_vec*np.log2(Pij_vec/(Pi_vec*Pj_vec)))
            mic_out["si"].append(i)
            mic_out["sj"].append(j)
            mic_out["mic"].append(mic)


    main_df = pd.DataFrame(main_out)
    mic_df = pd.DataFrame(mic_out)
    
    return main_df, mic_df