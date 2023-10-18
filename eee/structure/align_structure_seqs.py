
from eee.core.data import AA_3TO1

from eee._private.interface import create_new_dir
from eee._private.interface import launch
from eee._private.interface import rmtree
from eee._private.interface import WrappedFunctionException
from eee._private import logger

import numpy as np

import os

def _run_muscle(seq_list,
                muscle_binary="muscle",
                verbose=False,
                keep_temporary=False):
    """
    Actually run muscle.
    """

    tmp_dir = create_new_dir()
    input_fasta = "tmp-align_input.fasta"
    output_fasta = "tmp-align_output.fasta"

    # Write input fasta
    with open(os.path.join(tmp_dir,input_fasta),'w') as f:
        for i, s in enumerate(seq_list):
            f.write(f">seq{i}\n{''.join(s)}\n")

    # Try both versions of muscle command line (v. 5, v. 3)
    cmd_list = [[muscle_binary,"-align",input_fasta,"-output",output_fasta],
                [muscle_binary,"-in",input_fasta,"-out",output_fasta]]
    
    # Run muscle
    successful = False

    for cmd in cmd_list:

        try:
            launch(cmd=cmd,
                   run_directory=tmp_dir,
                   suppress_output=(not verbose))
        except (RuntimeError,FileNotFoundError,WrappedFunctionException):
            continue
    
        successful = True

    # If we did not actually do alignment, throw error
    if not successful:
        err = f"Alignment failed! Is a recent version of muscle in the path?\n"
        raise RuntimeError(err)

    # Read output fasta file
    output = {}
    with open(os.path.join(tmp_dir,output_fasta)) as f:
        
        for line in f:
            if line.startswith(">"):
                key = int(line[4:])
                output[key] = []
            else:
                output[key].extend(list(line.strip()))

    # Delete temporary files
    if not keep_temporary:
        rmtree(tmp_dir)

    # Get alignment columns for each site
    keys = list(output.keys())
    keys.sort()

    return [output[k] for k in keys]

def align_structure_seqs(dfs,
                         muscle_binary="muscle",
                         verbose=False,
                         keep_temporary=False):
    """
    Use muscle to align sequences from rcsb files, then do some clean up. 
    Renumber residues so they match between structures. If a site has a mixture
    of CYS and SER across structures, mutate the SER to CYS. 
    
    Parameters
    ----------
    dfs : list
        list of pandas dataframes containing structures
    muscle_binary : str, default="muscle"
        path to muscle binary
    verbose : bool, default=True
        whether or not to print out muscle output
    keep_temporary : bool, default=False
        do not delete temporary files
        
    Returns
    -------
    dfs : list
        list of pandas dataframes with structures updated to have the shared_fx,
        alignment_site, and identical_aa columns. 
    """

    seq_list = []
    for df in dfs:
        mask = np.logical_and(df.atom == "CA",df["class"] == "ATOM")
        this_df = df.loc[mask,:]        
        seq_list.append([AA_3TO1[aa] for aa in this_df["resid"]])

    output = _run_muscle(seq_list=seq_list,
                         muscle_binary=muscle_binary,
                         verbose=verbose,
                         keep_temporary=keep_temporary)
    

    # Convert output into column indexes and column contents. For example: 
    # MAST-
    # -ASTP
    # Would yield: 
    # + column_contents: [[0],[0,1],[0,1],[0,1],[1]]
    # + column_indexes: [[0,1,2,3],[1,2,3,4]]
    column_indexes = [[] for _ in range(len(output))]
    column_contents = [[] for _ in range(len(output[0]))]
    for i in range(len(output)):
        
        counter = 0            
        for j, c in enumerate(output[i]):
            if c != "-":
                column_contents[j].append(i)
                column_indexes[i].append(counter)
                counter += 1


    # Check sequence identity. identical_aa is True if the amino acids are the 
    # same at that site (or a mix of Ser and Cys), False if they differ. Gaps
    # do not count as different. 
    ser_to_cys = {}
    identical_aa = np.ones(len(column_contents),dtype=float)
    for i in range(len(column_contents)):

        struct_seen = column_contents[i]
        aa_seen = list(set([output[j][i] for j in struct_seen]))

        # All same
        if len(aa_seen) == 1:
            continue

        # If mixture of Cys and Ser seen, mutate SER --> CYS
        if set(aa_seen) == set(["C","S"]):
            ser_to_cys[i] = []
            for j in range(len(column_contents[i])):
                if output[j][i] == "S":
                    ser_to_cys[i].append(j)
        else:
            identical_aa[i] = False

    # Get lists of all CA atoms and residues
    residues = []
    for df in dfs:

        df["_resid_key"] = list(zip(df["chain"],df["resid"],df["resid_num"]))
        
        mask = np.logical_and(df.atom == "CA",df["class"] == "ATOM")
        this_df = df.loc[mask,:]

        residues.append(list(this_df["_resid_key"]))

    # Create an array indicating the fraction of structures sharing the site. 
    shared_column = np.zeros(len(column_contents),dtype=float)
    num_structures = len(dfs)
    for i in range(len(column_contents)):
        shared_column[i] = len(column_contents[i])/num_structures

    # Go through the alignment for each structure
    for i in range(len(column_indexes)):
        
        this_df = dfs[i]
        this_df["shared_fx"] = 0.0
        this_df["identical_aa"] = 0.0
        
        for j in range(len(column_indexes[i])):
            
            idx = column_indexes[i][j]

            # Record how many structures share the column
            this_resid = residues[i][j]
            this_resid_mask = this_df["_resid_key"] == this_resid

            # Change residue numbering
            this_df.loc[this_resid_mask,"resid_num"] = f"{idx + 1:d}"

            # Record shared fraction
            this_df.loc[this_resid_mask,"shared_fx"] = shared_column[idx]

            # Record identical amino acids. 
            this_df.loc[this_resid_mask,"identical_aa"] = identical_aa[idx]        

            # Mutate ser to cys from sites with mix of ser and cys across the
            # structures
            if j in ser_to_cys:
                if i in ser_to_cys[j]:
                    if verbose:
                        logger.log(f"Introducing S{j}C into structure {i}")
                        
                    this_df.loc[this_resid_mask,"resid"] = "CYS"
                    atom_mask = np.logical_and(this_resid_mask,
                                               this_df["atom"] == "OG")
                    this_df.loc[atom_mask,"atom"] = "SG"

                    
    # Remove "_resid_key" convenience column
    for i in range(len(dfs)):
        dfs[i] = dfs[i].drop(columns="_resid_key")
        
    return dfs