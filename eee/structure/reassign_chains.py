from Bio import Align
import pandas as pd

from eee.structure.manipulate_chains import get_chain_df
from eee.structure.manipulate_chains import get_unique_chains
from eee.structure.manipulate_chains import chain_reindex
from eee._private import logger

def reassign_chains(dfs:list, ensemble:str,write_pdb=False, save_text_file=True):

    
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
        if row.iloc[0] not in unique_chains1:
            chain1 = chain1.drop(index=idx)
            
    
    #make a dictionary for chains that track the files that have this chain 
    #the dictionary has the chain ID from the first file as its key and a list of tuples as its value
    #the tuples are the df then the chain that is aligned with all others in the list
    dfs_by_chain={}
    for idx, row in chain1.iterrows():
        dfs_by_chain[row.iloc[0]]=[]
        chain_ID1=row.iloc[0]
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
            if row.iloc[0] not in unique_chains2:
                chain2 = chain2.drop(index=idx)
                
        
        #identify any chains that are 95% similar to each other
        #establish a count for these chains and raise an error if there is more than one chain in one structure matching with another
        #Given df1 with chain A and df2, only one chain in df2 should match with df1's chain A
        for idx1,row1 in chain1.iterrows():
            chain_ID1=row1.iloc[0]
            count=0
            for idx2,row2 in chain2.iterrows():
                chain_ID2=row2.iloc[0]
                scores=[]

                score=aligner.score(chain1.seq[idx1],chain2.seq[idx2])
                match=score/min(len(chain1.seq[idx1]),len(chain2.seq[idx2]))
                #if the match is sufficient we append a tuple with the name of the file and the chain of the file
                #into the dfs_by_chain dict
                if match >= 0.95 and count<1:
                    dfs_by_chain[chain_ID1].append((df2,chain_ID2))
                    count+=1
                elif match >= 0.95 and count >= 1:
                    df_name=df2.name[0]
                    logger.log(df_name+' has chains that are <95 percent similar to each other but 95 percent similar to other chains in the ensemble.')

                    

    #creating truly ubiquitous chain argument and setting to false
    truly_ubiquitous_chain=False
    #creating a dictionary that only includes chains that are present in all files in the ensemble
    ubiquitous_chains={}
    number_dfs_by_chain=[]
    for item in dfs_by_chain:
        number_dfs_by_chain.append(len(dfs_by_chain.get(item)))
        if len(dfs_by_chain.get(item)) == len(dfs):
            ubiquitous_chains[item]=dfs_by_chain.get(item)
            truly_ubiquitous_chain=True
                
            
    #if there is no chain that is in every structure, I will find the chain(s) present in the most structures and keep those
    if truly_ubiquitous_chain==False:
        most_structs=max(number_dfs_by_chain)
        for item in dfs_by_chain:
            if len(dfs_by_chain.get(item)) == most_structs:
                ubiquitous_chains[item]=dfs_by_chain.get(item)

    #create alphabet string and set idx for that string to 0
    alphabet='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    idx=0

    dfs=[]
    chains_for_thermo=[]
    #for any chain that is in all pdbs, rename chain so they are consistent with each other
    #start from A and move forward from there
    for c in ubiquitous_chains:
        chain=alphabet[idx]
        chains_for_thermo.append(chain)
        for f in ubiquitous_chains.get(c):
            if f[1] != chain:
                #the chain reindexing will impact the dfs that are in ubiquitous chains
                #so I don't need to add it to a new list yet
                chain_reindex(df=f[0],prev_chain=f[1],new_chain=chain)
        idx=idx+1
        
    #now add all the structures in the first chain of ubiquitous chains into a list
    for i in ubiquitous_chains[list(ubiquitous_chains.keys())[0]]:
        dfs.append(i[0])
        
    if save_text_file==True:
        #create a text file with each chain on its own line
        chains='\n'.join(chains_for_thermo)
        f = open(ensemble+'_chains.txt', 'w')
        f.write(chains)
        f.close()

    return dfs