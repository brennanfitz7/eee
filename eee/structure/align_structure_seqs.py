from eee.core.data import AA_3TO1

from eee._private.interface import create_new_dir
from eee._private.interface import launch
from eee._private.interface import rmtree
from eee._private.interface import WrappedFunctionException
from eee._private import logger

import numpy as np
import pandas as pd
from Bio import Align

import os
import statistics

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


def align_structure_seqs(original_dfs,
                         muscle_binary="muscle",
                         verbose=False,
                         keep_temporary=False):
    """
    Use muscle to align sequences from rcsb files, then do some clean up. 
    Renumber residues so they match between structures. If a site has a mixture
    of CYS and SER across structures, mutate the SER to CYS. 
    
    Parameters
    ----------
    original_dfs : list
        list of pandas dataframes containing structures of pdbs that have been through sync_structures
        (with pdb name as a column in the df)
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
    
    #setting up pairwise aligner
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -0.5
    aligner.extend_gap_score = -0.1
    aligner.target_end_gap_score = 0.0
    aligner.query_end_gap_score = 0.0

    pdb_list=[]
    original_seq_list=[]
    
    for df in original_dfs:
        #add to pdb_list
        pdb_list.append(df.name[1])
        #create a mask to only have one residue
        mask = np.logical_and(df.atom == "CA",df["class"] == "ATOM")
        this_df = df.loc[mask,:]  
        original_seq_list.append(''.join([AA_3TO1[aa] for aa in this_df["resid"]]))


    #make sequence df
    seq_dict = {'pdb': pdb_list, 'seq': original_seq_list} 
    seq_df = pd.DataFrame(seq_dict)  
    
    dfs=[]
    
    for i in range(0,len(seq_df)):
        matches=[]
        #now we do a pairwise alignment of this seq with every other seq
        for n in range(0,len(seq_df)):
            if i==n:
                continue
            else:
                score=aligner.score(seq_df.seq[i],seq_df.seq[n])
                match=score/min(len(seq_df.seq[i]),len(seq_df.seq[n]))
                matches.append(match)
            
        #get the mean match
        mean_match=statistics.mean(matches)
        
        if mean_match <= 0.90:
            print(seq_df.pdb[i]+' had a mean match score of '+str(mean_match)+'. It is being discarded.')
            
        if mean_match > 0.90:
            df=original_dfs[i]
            if seq_df.pdb[i] in df.name[1]:
                dfs.append(df)
            else:
                print('The dataframes and sequences are not in the same order')
                
    seq_list=[]
    
    for df in dfs:
        #create a mask to only have one residue
        mask = np.logical_and(df.atom == "CA",df["class"] == "ATOM")
        this_df = df.loc[mask,:]  
        seq_list.append(''.join([AA_3TO1[aa] for aa in this_df["resid"]]))


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
            column_indexes[i].append(counter)
            counter += 1
            if c != "-":
                column_contents[j].append(i)

    # Check sequence identity. identical_aa is True if the amino acids are the 
    # same at that site (or a mix of Ser and Cys), False if they differ. Gaps
    # do not count as different. 
    ser_to_cys = {}
    identical_aa = np.ones(len(column_contents),dtype=float)
    shared_column = np.zeros(len(column_contents),dtype=float)
    
    #Use column_indexes and column_contents to create a list indicating where residues are identical
    #between all structures
    for i in range(len(column_contents)):

        struct_seen = column_contents[i]

        aa_seen = list(set([output[j][i] for j in struct_seen]))

        if len(aa_seen) == 1 and len(column_contents[i])==len(column_indexes):
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
    for i in range(len(output)):
        
        #establishing which df we are editing and then creating shared_fx and identical_aa column to be 0
        this_df = dfs[i]
        this_df["shared_fx"] = 0.0
        this_df["identical_aa"] = 0.0
        
        #establish your various counters:
        #C will be the resid_number for residues shared between all structures
        #N will be the resid_number for residues that are not share between all structures
        #idx is what you will use to iterate over the residue list
        C=1
        N=-1
        idx=0
        
        for j in range(len(output[i])):

            #if there is a gap in the alignment, this will move down the lists based on the alignment
            #but will not affect the residues or residue numbers
            if output[i][j]=='-':
                continue
            
            else:
                # Create the mask that will select all rows for this residue
                this_resid = residues[i][idx]
                this_resid_mask = this_df["_resid_key"] == this_resid
                idx=idx+1
                
                # Record identical amino acids. 
                this_df.loc[this_resid_mask,"identical_aa"] = identical_aa[j]  

                #change identical_aa to a list
                list_identical_aa=list(identical_aa)
                
                #if a residue is shared between all structures, number it then add 1
                if list_identical_aa[j]==1:
                    this_df.loc[this_resid_mask,"resid_num"] = C
                    C=C+1
                
                #if a residue is not shared between all structures, number it negatively then add
                elif list_identical_aa[j]==0:
                    this_df.loc[this_resid_mask,"resid_num"]= N
                    N=N-1
                    
                    
                # Record shared fraction
                this_df.loc[this_resid_mask,"shared_fx"] = shared_column[j]


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


    #Remove "_resid_key" convenience column
    for i in range(len(dfs)):
        dfs[i] = dfs[i].drop(columns="_resid_key")

    return dfs