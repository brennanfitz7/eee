import pandas as pd



def convert_to_df(ddg_output:str):

    #1abc_CALCULATOR_raw_ddgs.txt
    """
    Creates a dataframe with each mutation and the resulting DDG.

    Parameters
    ----------
    
    ddg_output : str
        output of raw ddg calculations

    Returns
    -------
    Pandas dataframe with each mutation and the resulting DDG's.
    """
    
    pdb_id=ddg_output.split('_')[0]
    muts_file=pdb_id+'_acdc_muts.tsv'
    output_df=pdb_id+'_acdc_ddg_df.csv'
    
    with open(ddg_output,'r') as ddgs:
        output_list = ddgs.read().splitlines()
    
    ddg_list=[]
    for item in output_list:
        if item[3:len(item)].isdigit()==True:
            #adds tod ddg_list and changes sign so that stabilizing=neg, destab=positive
            ddg_list.append(-(float(item)))
        else:
            continue
    
    
    input_df=pd.read_table(muts_file, delimiter='\t',names=['Mutation', 'Prof', 'PDB','Chain'])
    
    mut_and_chain=[]
    for index,row in input_df.iterrows():
        mutation=input_df['Mutation'][index]
        chain=input_df['Chain'][index]
        mut_and_chain.append(mutation[0]+chain+mutation[1:])
        
    
    ddg_df = pd.DataFrame(
    {'Mutation': mut_and_chain,'DDG': ddg_list})

    ddg_df.to_csv(output_df, index=False)
    
    return ddg_df





