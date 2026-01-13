import json
import glob
import statistics
from Bio import Align

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from eee.structure.align_structure_seqs import _run_muscle
from eee._private import logger


def filter_seq_length(seq_df):

    """
     Filters out the sequences in the ensemble that are not the same size as the majority of the sequences. 

    Parameters
    ----------
        
    seq_df : pandas DataFrame
        a dataframe with columns "pdb" and "seq" containing all the pdbs in the ensemble.
        
    Returns
    -------
    seq_df only with sequences of the selected length
    """
    
    #create a column in the df with lengths of sequences
    len_seqs=[len(i) for i in seq_df.seq]
    
    seq_df['seq_length']=len_seqs
    
    #remove any sequences under 50 residues long
    seq_df=seq_df.loc[seq_df.seq_length>=50]
    
    seq_df=seq_df.reset_index(drop=True)
    
    #if all sequences are within 50 residues of each other, keep everything
    if max(seq_df.seq_length)-min(seq_df.seq_length)<50 :
        
        logger.log('All sequences are withing 50 residues length of each other.')
        
    else:
        
        logger.log('Filtering sequences by length.')
        
        #use the histogram function (binned by 5) to find the frequency of seq lengths
        hist_output=plt.hist(seq_df.seq_length, bins=range(min(seq_df.seq_length),max(seq_df.seq_length),5))
        
        bin_start=hist_output[1][0:-1]
        
        hist_df=pd.DataFrame({'counts':list(hist_output[0]),'bin_start':bin_start})
        
        #find the highest peak
        max_idx = hist_df.loc[hist_df['counts'] == max(hist_df.counts)].index[0]
        
        to_collect=[]
        
        #if the max_idx is not in the middle of the graph, set the stop and start to be within index
        #if max_idx is in the middle, look at the five peaks on either side
        if max_idx-5<0:
            
            start=0
            
        else:
            start=max_idx-5
            
        if max_idx+5 > (hist_df.index.stop -1):
            
            stop = hist_df.index.stop -1
            
        else:
            stop=max_idx+5
        
        #find peaks within this range that are at least 20% as large as the max peak
        for i in range(start,stop):

            if hist_df.counts[i] >= ((max(hist_df.counts))*0.2):

                to_collect.append(hist_df.loc[hist_df['counts']==hist_df.counts[i]].index[0])
        
        start_len=float(hist_df.loc[to_collect[0],'bin_start'])
        
        end_len=float(hist_df.loc[to_collect[-1],'bin_start']+5)
        
        #apply a mask to only select for sequences within the desired range
        mask=np.logical_and(seq_df.seq_length>=start_len,seq_df.seq_length<=end_len)
        
        seq_df=seq_df.loc[mask,:]
        
        seq_df=seq_df.reset_index(drop=True)
        
    seq_df=seq_df.drop('seq_length',axis=1)
        
    return seq_df



def align_thermo_output_seqs(original_dfs, 
                             muscle_binary="muscle",
                             verbose=False,
                             keep_temporary=False):
    #RUN THIS AFTER SEQ FILES ARE CONCATENATED SO DDG CSVS HAVE THE SAME FILE
    #RUN BEFORE REMOVAL OF SELF TO SELF MUTATIONS
    
    #setting up pairwise aligner
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -0.5
    aligner.extend_gap_score = -0.1
    aligner.target_end_gap_score = 0.0
    aligner.query_end_gap_score = 0.0

    pdb_list=[]
    seq_list=[]
    
    for df in original_dfs:
        #add to pdb_list
        pdb_list.append(df.pdb[1].split('_')[0])
        #create a mask to only have one residue
        mask = np.array(df.mutation=='A')
        this_df = df.loc[mask,:]  
        seq_list.append(''.join(this_df["wildtype"]))

    #make sequence df
    seq_dict = {'pdb': pdb_list, 'seq': seq_list} 
    seq_df = pd.DataFrame(seq_dict) 

    #filter by seq_length to make sure all sequences are approx the same size
    seq_df=filter_seq_length(seq_df)
    
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
            if seq_df.pdb[i] in df.pdb[1]:
                dfs.append(df)
            else:
                print('The dataframes and sequences are not in the same order')
                
    seq_list=[]
    
    for df in dfs:
        #create a mask to only have one residue
        mask = np.array(df.mutation=='A')
        this_df = df.loc[mask,:]  
        seq_list.append(this_df["wildtype"])
        
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
    # same at that site , False if they differ. Gaps do not count as different. 
    identical_aa = np.ones(len(column_contents),dtype=float)
    shared_column = np.zeros(len(column_contents),dtype=float)
    
    #Use column_indexes and column_contents to create a list indicating where residues are identical
    #between all structures
    for i in range(len(column_contents)):

        struct_seen = column_contents[i]

        aa_seen = list(set([output[j][i] for j in struct_seen]))

        if len(aa_seen) == 1 and len(column_contents[i])==len(column_indexes):
            continue
            
        else:
            identical_aa[i] = False
            
            
    # Get list of residues by selecting only one mutation
    residues = []
    for df in dfs:

        df["_resid_key"] = list(zip(df["chain"],df["wildtype"],df["position"]))
        
        mask = np.array(df.mutation=='A')
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
                    this_df.loc[this_resid_mask,"position"] = C
                    C=C+1
                
                #if a residue is not shared between all structures, number it negatively then add
                elif list_identical_aa[j]==0:
                    this_df.loc[this_resid_mask,"position"]= N
                    N=N-1
                    
                    
                # Record shared fraction
                this_df.loc[this_resid_mask,"shared_fx"] = shared_column[j]



    #Remove "_resid_key" convenience column
    for i in range(len(dfs)):
        dfs[i] = dfs[i].drop(columns="_resid_key")

    return dfs




def get_combined_df(folder:str,prot_name:str,need_name_dict:bool,unnamed_col_exists:bool, use_ddg_mult:bool):
    #FIX DDG MULT BEFORE I USE IT AGAIN

    """
     Compiles ensemble results from ThermoMPNN. 

    Parameters
    ----------
        
    folder : str
        folder that contains ThermoMPNN ddg results, multiplier dictionaries, and a name dictionary for all structures in an ensemble.

    prot_name : str
        name of the protein

    need_name_dict : bool
        specifies whether these proteins need a name dict or if they will just be categorized under their pdb file name
        
    Returns
    -------
    the ensemble df
    """
    
    if need_name_dict==True:
        #getting name_dict from json
        with open(folder+'/name_dict.json', 'r') as openfile:
                    name_dict = json.load(openfile)
            
    #getting list of ddgs
    ddg_list=glob.glob(folder+'/*.csv')

    #making lists
    single_chain_dfs=[]
    name_list=[]
    df_list=[]
    
    #creating "raw_dfs" for each of the ddg_csvs and assigning them names
    for item in ddg_list:
        pdb_id=item.split('/')[-1].split('_')[0]
        single_df = pd.read_csv(item)
        single_chain_dfs.append(single_df)
        single_df.name=pdb_id
        
        #adding pdb_ids to the name_list
        if single_df.name not in name_list:
            name_list.append(single_df.name)
            
    #make list of dfs for input into align_thermo_output_seqs
    raw_df_list=[]
            
    #sorting dfs by the pdb_id (now known as tags) 
    for tag in name_list:
        same_pdb_list=[]
        for df in single_chain_dfs:
            my_tag=df.name
            if my_tag==tag:
                same_pdb_list.append(df)
        #concatenating all dfs from the same pdb then sorting them 
        raw_df=pd.concat(same_pdb_list)
        raw_df.sort_values(by=['chain','position'],inplace=True)
        raw_df.reset_index(inplace=True)
        raw_df.name=tag
        raw_df_list.append(raw_df)
        
    
    #run align thermo output seqs to get the numbering right
    aligned_dfs=align_thermo_output_seqs(raw_df_list)
    
    #set up df list to put new_dfs in
    df_list=[]
    
    for aligned_df in aligned_dfs:
        
        #get the pdb id before deleting column
        name=aligned_df.pdb[0].split('.')[0]

        #drop any self to self mutations
        for idx,row in aligned_df.iterrows():
            if aligned_df.wildtype[idx]==aligned_df.mutation[idx]:
                aligned_df.drop([idx],inplace=True)
        
        #drop any residues with a negative number in position (not shared by other structures)
        for idx,row in aligned_df.iterrows():
            if aligned_df.position[idx] <= 0:
                aligned_df.drop([idx],inplace=True) 
        
        if use_ddg_mult==True:
            #use ddg_mults to multiply ddg values
            #THIS CODE NEEDS TO BE INVESTIGATED AND FIXED BEFORE I USE IT AGAIN
            with open(folder+'/'+pdb_id+'_ddg_mult.json', 'r') as openfile:
                mult_dict = json.load(openfile)
            for idx,row in aligned_df.iterrows():
                res_chain=aligned_df.chain[idx]
                if aligned_df.chain[idx] in mult_dict.keys():
                    aligned_df.loc[idx, "ddG_pred"] = mult_dict[res_chain]*aligned_df[idx,'ddG_pred']
                
        #create column with wildytpe, residue number, mutation
        mut_col=[]
        for idx,row in aligned_df.iterrows():
            mut=aligned_df.wildtype[idx]+str(aligned_df.position[idx])+aligned_df.mutation[idx]
            mut_col.append(mut)
        
        aligned_df['mut']=mut_col
        
        #drop columns that are no longer needed 
        if unnamed_col_exists==True:
            new_df=aligned_df.drop(axis=1,labels=['Unnamed: 0','Model','Dataset','pdb','shared_fx','identical_aa','index','mutation','wildtype','chain'])
        elif unnamed_col_exists==False:
             new_df=aligned_df.drop(axis=1,labels=['Model','Dataset','pdb','shared_fx','identical_aa','index','mutation','wildtype','chain'])

               
        #name the df
        new_df.name=name

        
        #append to df_list
        df_list.append(new_df)

    #now time to combine the dfs
    mut_sets=[]
    for df in df_list:
        mut_sets.append(set(df.mut))

        #create a set of all shared mutations
        shared_muts=mut_sets[0].intersection(*mut_sets[1:])
        
    
    muts=list(shared_muts)
    
    combined_df = pd.DataFrame(muts, columns=['mut'])

    for df in df_list:
        name=df.name
        
        #creating placeholder column full of ones with the correct name
        placeholder = np.ones(len(muts),dtype=float)
        combined_df[name]=placeholder
        
        #now replace these ones with the ddG_pred values
        
        for idx,row in combined_df.iterrows():
            
            mut_site = combined_df.mut[idx]

            ddG_series = df.loc[df.mut==mut_site,'ddG_pred']
            
            mut_ddG=ddG_series.iloc[0]

            combined_df.loc[idx,name] = mut_ddG
            
                 
    #adding site and unfolded columns to combined_df
    combined_df['unfolded']=np.zeros(len(combined_df),dtype=float)
    
    site=[]
    for item in combined_df.mut:
        site.append(float(item[1:-1]))

    combined_df.insert(loc=0, column='site', value=site)
    
    #sort combined_df to get it all in order
    combined_df.sort_values(by=['site','mut'],inplace=True)
    
    combined_df.to_csv(folder+'/'+prot_name+'_thermo_ddgs.csv',index=False)

    return combined_df
    
