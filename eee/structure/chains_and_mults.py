from eee.io.read_structure import read_structure
from eee.core.data import AA_3TO1

import pandas as pd
import numpy as np
import json

from Bio import Align

def chains_and_multipliers(pdb_file, save_dict:bool):
    """
    Find unique chains and multipliers within a single pdb file
        
    Parameters
        ----------
        pdb_file : str
            pdb for which you are finding chains and multipliers    
        save_dict : bool
            bool about whether to save the mult_dict as json file or not

        Returns
        -------
        mult_dict : dictionary
            a dictionary where the keys are unique chains within the pdb_file and the values are the number of copies of this chain in the pdb
        """ 
    #establish pdb_id
    pdb_id=pdb_file.split('.')[0]
    
    #get dataframe from pdb, select for atoms, and drop duplicates
    raw_prot=read_structure(pdb_file)
    redundant_prot_seq=raw_prot[raw_prot['class'] == "ATOM"].loc[:,["chain",'resid','resid_num']].drop_duplicates(subset=["resid_num","chain"], keep='first')
    
    #getting rid of redundant chains
    #making a dataframe of Chain ID and sequence
    chain_df=pd.DataFrame(columns=['Chain_ID','seq'])
    for item in redundant_prot_seq.chain.unique():
        chain=str(item)
        temp_list=[]
        for index,row in redundant_prot_seq.iterrows():
            if chain==row[0]:
                temp_list.append(AA_3TO1.get(row[1]))

        test_string=''.join(temp_list)
        chain_df.loc[len(chain_df.index)]=[chain,test_string]
        
    #finding original chains 
    #setting up pairwise aligner
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -0.5
    aligner.extend_gap_score = -0.1
    aligner.target_end_gap_score = 0.0
    aligner.query_end_gap_score = 0.0

    #now time to do alignment between each different chain
    redundant=[]
    mult_dict={}
    for i in range(0,len(chain_df)):
        #this line is to make sure we don't run chains that have already been grouped with other chains
        if chain_df.Chain_ID[i] not in redundant:
            main_chain=chain_df.Chain_ID[i]
            #creating list with this chain ID and all chains like it
            similar_seqs=[main_chain]
            #now we do a pairwise alignment of this chain with every other chain (that hasn't been grouped already)
            for n in range(i+1,len(chain_df)):
                scores=[]
                for alignment in aligner.align(chain_df.seq[i],chain_df.seq[n]):
                    scores.append(alignment.score)
                mean_score=np.mean(scores)
                match=mean_score/min(len(chain_df.seq[i]),len(chain_df.seq[n]))
                #if the match is sufficient we append the chain to the similar seqs list
                if match >= 0.99:
                    redundant.append(chain_df.Chain_ID[n])
                    similar_seqs.append(chain_df.Chain_ID[n])
            mult_dict[main_chain]=len(similar_seqs)

    if save_dict==True:        
        #now write out the dict to the json
        with open(pdb_id+"_ddg_mult.json", "w") as outfile:
            json.dump(mult_dict, outfile)
    
    return mult_dict