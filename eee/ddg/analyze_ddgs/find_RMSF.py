from eee.io.read_structure import read_structure
from eee._private import logger
from eee.structure.align_structure_seqs import align_structure_seqs
from eee.structure.align_structures import align_structures

import pandas as pd
import numpy as np

import statistics
import glob


def get_RMSF(x_values,y_values,z_values):
    
    #check that all value lists are the same length
    if len(x_values) == len(y_values) == len(z_values):
        pass
    else:
        err = "Value lists must all be the same length!"
        raise ValueError(err)
    
    x_mean=statistics.mean(x_values)
    y_mean=statistics.mean(y_values)
    z_mean=statistics.mean(z_values)
    
    squared_diff=[]
    
    for i in range(0,len(x_values)):
        squared_diff.append(((x_values[i]-x_mean)**2)+((y_values[i]-y_mean)**2)+((z_values[i]-z_mean)**2))
    
    RMSF=(((sum(squared_diff))/len(x_values))**0.5)
    
    return RMSF
        

def get_RMSF_ens(ens_dir, chains_file):
    
    """
    ens_dir : str
        this ensemble directory should contain pdbs that have been through sync_structures
    
    chains_file : str
        file containing the chains that are in all pdbs in the ensemble
    """
    
    pdb_list=glob.glob(ens_dir+'/*pdb')
    
    #read synced_chains file
    chains=[]
    with open(chains_file, 'r') as file:
        for line in file:
            chains.append(line.strip())
    
    dfs=[]
    
    #read the pdb structures and add a name column
    logger.log('Reading in pdbs.')
    for item in pdb_list:
        try:
            read_df=read_structure(item,remove_multiple_models=False)
            read_df['name']=item
            dfs.append(read_df)
        except:
            logger.log('PDB '+str(item)+' failed when reading structure and incorporating models.')

    #align the structure sequences
    logger.log('Aligning structure sequences.')
    aligned_dfs=align_structure_seqs(original_dfs=dfs,limited_chains=chains)
    
    #align the dfs by chain
    #if there are multiple chains in synced_chains, find the longest one
    if len(chains)>1:
        chain_lengths={}
        for chain in chains:
            length=[]
            for df in aligned_dfs:
                length.append(len(df.loc[df['chain']==chain]))
                
            chain_lengths[np.mean(length)]=chain
        
        mychain=chain_lengths.get(max(chain_lengths.keys()))
    
    else:
        mychain=chains[0]
        
    
    #align structures with lovoalign
    logger.log('Aligning structures by chain with lovoalign.')
    aligned_dfs = align_structures(aligned_dfs,chain=mychain)
    
            
    prepped_dfs=[]
    resid_sets=[]
    
    #prep the dfs by creating a resid key column and keeping only alpha carbons and regular atoms
    for df in aligned_dfs:

        df["_resid_key"] = list(zip(df["chain"],df["resid"],df["resid_num"]))

        mask=np.logical_and(df['atom']=='CA',df['class']=='ATOM')
        this_df = df.loc[mask,:]

        resid_sets.append(set(this_df._resid_key))
        
        prepped_dfs.append(this_df)

    #create a list of all shared residues
    shared_resids=list(resid_sets[0].intersection(*resid_sets[1:]))

    #make sure that all of the dfs have only shared residues
    only_shared=[]
    for df in prepped_dfs:
        this_df=df.loc[df['_resid_key'].isin(shared_resids)]
        only_shared.append(this_df)
        
    RMSF_df = only_shared[0][['chain','resid','resid_num','_resid_key']].copy()
    
    #go through RMSF df, collect values for each site from each dataframe, and get the RMSF value
    logger.log('Finding RMSF values for each site.')
    RMSF=[]
    for idx, row in RMSF_df.iterrows():

        resid_key=RMSF_df.loc[idx,'_resid_key']

        x_vals=[]
        y_vals=[]
        z_vals=[]

        #get all my x, y, and z values in my file
        for df in only_shared:
            x_vals.append(float(df.loc[df['_resid_key']==resid_key,'x']))
            y_vals.append(float(df.loc[df['_resid_key']==resid_key,'y']))
            z_vals.append(float(df.loc[df['_resid_key']==resid_key,'z']))


        RMSF.append(get_RMSF(x_vals, y_vals, z_vals))
        
    #add new RMSF column to my df
    RMSF_df.insert(3, "RMSF", RMSF)
    
    #save df as a csv
    dir_name=ens_dir.split('/')[-1]
    ens_alone=dir_name.split('_synced')[0]
    RMSF_df.to_csv(ens_dir+'/'+ens_alone+'_RMSF.csv',index=False)
    
    return RMSF_df