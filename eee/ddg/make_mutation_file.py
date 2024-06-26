from eee.io.read_structure import read_structure

import pandas as pd
import json

def make_mutation_file(pdb_file:str):
    """
    Makes a dataframe with all possible mutations from a pdb file. 

    Parameters
    ----------
        
    pdb_file : str
        pdb file of protein of interest

    Returns
    -------
    Dataframe with all possible mutations from a pdb file. 
    """

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
    original_chains=[]
    all_chains=[]
    for group, df in chain_df.groupby('seq'):
        temp_list=list(df['Chain ID'])
        original_chains.append(temp_list[0])
        all_chains.append(temp_list)

    #creating ddg multiplier json
    mult_dict={}
    for item in all_chains:
        mult_dict[item[0]]=len(item)

    with open(pdb_id+"_ddg_mult.json", "w") as outfile:
        json.dump(mult_dict, outfile)
        
    prot_seq=redundant_prot_seq[(redundant_prot_seq['chain'].isin(original_chains))]

    #create new mutation dictionary
    mutation_dict={'chain':[],'wt_res':[],'res_pos':[],'mut_res':[]}
    
    aa_list=["ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR"]
    
    #populate mutation dictionary with info from original df
    for index, row in prot_seq.iterrows():
        for aa in aa_list:
            if row[1]!=aa:
                mutation_dict["chain"].append(row[0])
                mutation_dict["wt_res"].append(row[1])
                mutation_dict["res_pos"].append(row[2])
                mutation_dict["mut_res"].append(aa)
            else:
                continue
    
    #create a dataframe from the mutation dictionary
    mutation_df=pd.DataFrame.from_dict(mutation_dict)
    
    return mutation_df