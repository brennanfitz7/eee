from eee.io.read_structure import read_structure

import pandas as pd
import json

def chains_and_multipliers(pdb_file):
    #establish pdb_id
    pdb_id=pdb_file.split('.')[0]
    
    #get dataframe from pdb, select for atoms, and drop duplicates
    raw_prot=read_structure(pdb_file)
    redundant_prot_seq=raw_prot[raw_prot['class'] == "ATOM"].loc[:,["chain",'resid','resid_num']].drop_duplicates(subset=["resid_num","chain"], keep='first')
    
    #getting rid of redundant chains
    #making a dataframe of Chain ID and sequence
    chain_df=pd.DataFrame(columns=['Chain ID','seq'])
    for item in redundant_prot_seq.chain.unique():
        chain=str(item)
        temp_list=[]
        for index,row in redundant_prot_seq.iterrows():
            if chain==row[0]:
                temp_list.append(row[1])

        test_string=','.join(temp_list)
        chain_df.loc[len(chain_df.index)]=[chain,test_string]
        
    #finding original chains 
    all_chains=[]
    for group, df in chain_df.groupby('seq'):
        temp_list=list(df['Chain ID'])
        all_chains.append(temp_list)

    #finding original chains 
    all_chains=[]
    for group, df in chain_df.groupby('seq'):
        temp_list=list(df['Chain ID'])
        all_chains.append(temp_list)

    #creating ddg multiplier json
    mult_dict={}
    for item in all_chains:
        mult_dict[item[0]]=len(item)
    
    with open(pdb_id+"_ddg_mult.json", "w") as outfile:
        json.dump(mult_dict, outfile)
    
    return mult_dict