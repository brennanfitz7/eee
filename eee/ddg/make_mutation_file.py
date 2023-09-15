from eee.io import read_structure

import pandas as pd

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
    
    #get dataframe from pdb, select for atoms, and drop duplicates
    raw_prot=read_structure(pdb_file)
    prot_seq=raw_prot[raw_prot['class'] == "ATOM"].loc[:,["chain",'resid','resid_num']].drop_duplicates(subset=["resid_num","chain"], keep='first')
    
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