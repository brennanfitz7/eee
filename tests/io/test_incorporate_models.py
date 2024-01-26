from eee.io.read_structure import read_structure

import pandas as pd


def test_run_ensemble(test_pdbs):
    pdb_dfs=[]
    for pdb in test_pdbs:
        pdb_dfs.append(read_structure(pdb, remove_multiple_models=False))

    #test that alphabet list is created correctly
    alphabet_list=[x.upper() for x in (list(map(chr, range(ord('a'), ord('z')+1))))]
    assert len(alphabet_list)==26
    for item in alphabet_list:
        assert item.isupper()

    for df in pdb_dfs:
        
        dfs_by_model=[]
        assert issubclass(type(dfs_by_model), list)

        for model in df['model'].unique():
            grouped = df.groupby(df.model)
            df_new = grouped.get_group(model)
            dfs_by_model.append(df_new)

        #test that the number of dfs you end up with is the number of models in original pdb
        assert len(dfs_by_model) == len(df['model'].unique())

        df_new_chains=[]
        assert issubclass(type(df_new_chains),list)

        for df in dfs_by_model:
            #create list of unique chain ID's in a df
            old_chains=list(df['chain'].unique())
            
            for item in old_chains:
                #remove any letter used in chain so it can't be reused
                if item in alphabet_list:
                    alphabet_list.remove(item)
                    #test that item is removed
                    assert item not in alphabet_list
                #rename chain if it has already been used
                elif item not in alphabet_list:
                    new_chain_ID=alphabet_list[0]
                    df.loc[df['chain'] == item, 'chain'] = new_chain_ID
                    #assert that new_chain_ID is in the df and the old one is not
                    assert new_chain_ID in list(df['chain'])
                    assert item not in list(df['chain'])
                    #remove this new chain ID so it can't be reused
                    alphabet_list.remove(new_chain_ID)
                    #test that item is removed
                    assert new_chain_ID not in alphabet_list

                        #populate list of dfs with new chains
        
            df_new_chains.append(df)
        
        #test that the number of new dfs is consistent with the number of models
        assert len(df['model'].unique()) ==len(df_new_chains)

        #concatenate all dfs with new chains into one df
        df_unique_chains=pd.concat(df_new_chains)

        #check that the original df and old df are the same length
        assert len(df) == len(df_unique_chains)
        assert issubclass(type(df_unique_chains), pd.DataFrame)
