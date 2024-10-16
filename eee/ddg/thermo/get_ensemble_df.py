import pandas as pd

import json
import glob

def get_ensemble_df(folder:str,prot_name:str,need_name_dict:bool,unnamed_col_exists:bool):

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
            
    #sorting dfs by the pdb_id (now known as tags) 
    for tag in name_list:
        same_pdb_list=[]
        for df in single_chain_dfs:
            my_tag=df.name
            if my_tag==tag:
                same_pdb_list.append(df)
        #concatenating all dfs from the same pdb then sorting them 
        raw_df=pd.concat(same_pdb_list)
        raw_df.sort_values(by=['pos'],inplace=True)
        raw_df.reset_index(inplace=True)
        raw_df.name=tag
        
        #drop any self to self mutations
        for idx,row in raw_df.iterrows():
            if raw_df.wtAA[idx]==raw_df.mutAA[idx]:
                raw_df.drop([idx],inplace=True)
        if unnamed_col_exists==True:
            new_df=raw_df.drop(axis=1,labels=['Unnamed: 0','pos','wtAA','mutAA'])
        elif unnamed_col_exists==False:
             new_df=raw_df.drop(axis=1,labels=['pos','wtAA','mutAA'])
        with open(folder+'/'+pdb_id+'_ddg_mult.json', 'r') as openfile:
                mult_dict = json.load(openfile)
        chain=item[-5]
        for index,row in new_df.iterrows():
            #my_mut=new_df.Mutation[index]
            #new_df.loc[index,'Mutation'] = my_mut[0]+chain+my_mut[1:]
            if chain in mult_dict.keys():
                new_df.loc[index, "ddG (kcal/mol)"] = mult_dict[chain]*new_df['ddG (kcal/mol)'][index]
                new_df.name=tag
        df_list.append(new_df)
        
    mut_sets=[]
    for df in df_list:
        mut_sets.append(set(df.Mutation))

        #create a set of all shared mutations
        shared_muts=mut_sets[0].intersection(*mut_sets[1:])

    clean_df_list=[]
    for df in df_list:
        #make a mask to select for mutations shared between all files
        mask=df.Mutation.isin(shared_muts)
        #apply mask to df
        clean_df=df.loc[mask,:]
        #make sure the name of the df is maintained
        clean_df.name=df.name
        clean_df_list.append(clean_df)

    #making combined_df--starting with "mut" column
    mut_source=clean_df_list[0]
    mut_list=mut_source['Mutation'].tolist()
    combined_df=pd.DataFrame(mut_list, columns=['mut'])

    for df in clean_df_list:
        DDG_col = df['ddG (kcal/mol)'].tolist()
        if need_name_dict==True:
            #use name_dict to get the correct name for the pdb file
            name_of_df=name_dict.get(df.name)
        if need_name_dict==False:
             #use pdb file as name of column
             name_of_df=df.name
        combined_df[name_of_df] = DDG_col

    #making site column
    site=[]
    for index,row in combined_df.iterrows():
        site.append(row[0][1:-1])
    combined_df.insert(loc=0, column='site',value=site)

    #making unfolded column
    end=len(combined_df.columns)
    combined_df.insert(loc=end, column='unfolded',value=0)
    
    combined_df.to_csv(folder+'/'+prot_name+'_thermo_ddgs.csv', index=False)
    
    return combined_df