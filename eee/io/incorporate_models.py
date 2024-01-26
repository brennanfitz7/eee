"""
A function that changes multiple models into one model with many chains.
"""

import pandas as pd

def incorporate_models(df):
    """
    Take a pdb with multiple models that are used to represent
    a biological assembly and rename chains so that the models 
    are now obsolete and do not need to be maintained throughout 
    EEE pipeline.

    Parameters
    ----------
    input : pandas.DataFrame
       A pandas dataframe with information from a PDB file (obtained from read_structure.py)

    Returns
    -------
    pandas.DataFrame
        dataframe with chains renamed
    """

    #create alphabet list
    alphabet_list=[x.upper() for x in (list(map(chr, range(ord('a'), ord('z')+1))))]

    #split pdb df into separate dfs based on model
    dfs_by_model=[]

    for model in df['model'].unique():
        grouped = df.groupby(df.model)
        df_new = grouped.get_group(model)
        dfs_by_model.append(df_new)

    df_new_chains=[]

    for df in dfs_by_model:
        #create list of unique chain ID's in a df
        old_chains=list(df['chain'].unique())
        
        for item in old_chains:
            #remove any letter used in chain so it can't be reused
            if item in alphabet_list:
                alphabet_list.remove(item)
            #rename chain if it has already been used
            elif item not in alphabet_list:
                new_chain_ID=alphabet_list[0]
                df.loc[df['chain'] == item, 'chain'] = new_chain_ID
                #remove this new chain ID so it can't be reused
                alphabet_list.remove(new_chain_ID)
        
        #populate list of dfs with new chains
        df_new_chains.append(df)

    #concatenate all dfs with new chains into one df
    df_unique_chains=pd.concat(df_new_chains)

    return df_unique_chains