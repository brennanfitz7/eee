import pandas as pd



def convert_to_df(pdb_file:str):
    """
    Creates a dataframe with each mutation and the resulting DDG.

    Parameters
    ----------
    
    pdb_file : str
        pdb file of protein of interest

    Returns
    -------
    Pandas dataframe with each mutation and the resulting DDG's.
    """
    
    pdb_id=pdb_file.split('.')[0]
    output_file=pdb_id+'_ddg_output.txt'
    tsv_file=pdb_id+'_ddg_input.tsv'
    output_df=pdb_id+'_ddg_df.csv'
    
    with open(output_file,'r') as ddgs:
        output_list = ddgs.read().splitlines()
    
    ddg_list=[]
    for item in output_list:
        if "==" in item or len(item)<=3:
            continue
        else:
            ddg_list.append(item)
    
    input_df=pd.read_table(tsv_file, delimiter='\t',names=['Mutation', 'Prof', 'PDB','Chain'])
    
    ddg_df = pd.DataFrame(
    {'Mutation': input_df['Mutation'],'DDG': ddg_list})

    ddg_df.to_csv(output_df, index=False)
    
    return ddg_df



