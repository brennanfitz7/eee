import json
import glob
from Bio import Align
import pandas as pd
import numpy as np

from eee.io.read_structure import read_structure
from eee.io.write_pdb import write_pdb
from eee.core.data import AA_3TO1
from eee._private.logger import log
from eee.structure.chains_and_mults import chains_and_multipliers

def chain_reindex(df,prev_chain,new_chain,pdb=None,write_pdb=False):
    
    #first change chain names
    
    #if there is a previous chain with the new chain name I want, change it to a placeholder chain
    df.loc[df.chain == new_chain, 'chain'] = 'ph'
    
    #now change my prev_chain residues to new_chain
    df.loc[df.chain == prev_chain,'chain'] = new_chain
    
    #change the placeholder chain to the prev_chain (keeping it simple)
    df.loc[df.chain == 'ph','chain'] = prev_chain
    
    #then sort by chain (sort by atom number secondarily so as to not mess with the order of the atoms)
    
    df=df.sort_values(['chain', 'atom_num'], ascending=[True, True])
    
    if write_pdb==True:
        write_pdb(df, pdb, overwrite=True,write_with_models=True)
    
    return df

def get_chain_df(df):
    #select for atoms and drop duplicates
    redundant_prot_seq=df[df['class'] == "ATOM"].loc[:,["chain",'resid','resid_num']].drop_duplicates(subset=["resid_num","chain"], keep='first')

    #getting rid of redundant chains
    #making a dataframe of Chain ID and sequence
    chain_df=pd.DataFrame(columns=['Chain_ID','seq'])
    for item in redundant_prot_seq.chain.unique():
        chain=str(item)
        temp_list=[]
        for index,row in redundant_prot_seq.iterrows():
            if chain==row[0]:
                temp_list.append((AA_3TO1.get(row[1])))

        test_string=''.join(temp_list)
        chain_df.loc[len(chain_df.index)]=[chain,test_string]

    return chain_df


def get_unique_chains(df):

    #get df made up of chains
    chain_df=get_chain_df(df)
        
    #finding original chains 
    #setting up pairwise aligner
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -0.5
    aligner.extend_gap_score = -0.1
    aligner.target_end_gap_score = 0.0
    aligner.query_end_gap_score = 0.0

    #now time to do alignment between each different chain
    redundant=[]
    unique_chains=[]
    for i in range(0,len(chain_df)):
        #this line is to make sure we don't run chains that have already been grouped with other chains
        if chain_df.Chain_ID[i] not in redundant:
            main_chain=chain_df.Chain_ID[i]
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
            unique_chains.append(main_chain)


    return unique_chains

def reassign_chains(dfs:list, ensemble:str,write_pdb=False):
    
    #setting up pairwise aligner
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -0.5
    aligner.extend_gap_score = -0.1
    aligner.target_end_gap_score = 0.0
    aligner.query_end_gap_score = 0.0
    
    #starting with file 1--all other files will be compared to this one
    df1=dfs[0]
    
    #get the df sorted by chains for this first df
    chain1=get_chain_df(df1)
    
    #find unique chains and drop any chains that are redundant
    unique_chains1=get_unique_chains(df1)
    for idx, row in chain1.iterrows():
        if row[0] not in unique_chains1:
            chain1 = chain1.drop(index=idx)
            
            
    #make a dictionary for chains that track the files that have this chain 
    #the dictionary has the chain ID from the first file as its key and a list of tuples as its value
    #the tuples are the df then the chain that is aligned with all others in the list
    dfs_by_chain={}
    for idx, row in chain1.iterrows():
        dfs_by_chain[row[0]]=[]
        chain_ID1=row[0]
        dfs_by_chain[chain_ID1].append((df1,chain_ID1))
        
        
    #iterating over all other files
    for j in range(1,len(dfs)):
        
        #assign df2
        df2=dfs[j]
        
        #get the df sorted by chains for this first df
        chain2=get_chain_df(df2)
        
        #find unique chains and drop any chains that are redundant
        unique_chains2=get_unique_chains(df2)
        for idx, row in chain2.iterrows():
            if row[0] not in unique_chains2:
                chain2 = chain2.drop(index=idx)

        #identify any chains that are 95% similar to each other
        for idx1,row1 in chain1.iterrows():
            chain_ID1=row1[0]

            for idx2,row2 in chain2.iterrows():
                chain_ID2=row2[0]
                scores=[]

                for alignment in aligner.align(chain1.seq[idx1],chain2.seq[idx2]):
                    scores.append(alignment.score)
                mean_score=np.mean(scores)
                match=mean_score/min(len(chain1.seq[idx1]),len(chain2.seq[idx2]))
                #if the match is sufficient we append a tuple with the name of the file and the chain of the file
                #into the dfs_by_chain dict
                if match >= 0.95:
                    dfs_by_chain[chain_ID1].append((df2,chain_ID2))
                    
    #creating a dictionary that only includes chains that are present in all files in the ensemble
    ubiquitous_chains={}
    for item in dfs_by_chain:
        if len(dfs_by_chain.get(item)) == len(dfs):
            ubiquitous_chains[item]=dfs_by_chain.get(item)


    #create alphabet string and set idx for that string to 0
    alphabet='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    idx=0

    dfs=[]
    #for any chain that is in all pdbs, rename chain so they are consistent with each other
    #start from A and move forward from there
    for c in ubiquitous_chains:
        chain=alphabet[idx]
        for f in ubiquitous_chains.get(c):
            if f[1]==chain:
                continue
                dfs.append(f[0])
            else:
                df=chain_reindex(df=f[0],prev_chain=f[1],new_chain=chain)
                dfs.append(df)

        idx=idx+1
        
    #create a text file with each chain on its own line
    ubiq_list=list(ubiquitous_chains.keys())
    chains='\n'.join(ubiq_list)
    f = open(ensemble+'_chains.txt', 'w')
    f.write(chains)
    f.close()

    return dfs
