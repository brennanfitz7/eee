import pandas as pd

import json
import glob

def get_ensemble_df(folder:str,prot_name:str):

    """
     Compiles ensemble results from ThermoMPNN. 

    Parameters
    ----------
        
    folder : str
        folder that contains ThermoMPNN ddg results, multiplier dictionaries, and a name dictionary for all structures in an ensemble.

    prot_name : str
        name of the protein
        
    Returns
    -------
    the ensemble df
    """
    
    #getting name_dict from json
    with open(folder+'/name_dict.json', 'r') as openfile:
                name_dict = json.load(openfile)
            
    #getting list of ddgs
    ddg_list=glob.glob(folder+'/*.csv')
    
    df_list=[]
    for item in ddg_list:
        pdb_id=item.split('/')[-1].split('_')[0]
        raw_df = pd.read_csv(item)
        new_df=raw_df.drop(axis=1,labels=['Unnamed: 0','pos','wtAA','mutAA'])
        with open(folder+'/'+pdb_id+'_ddg_mult.json', 'r') as openfile:
                mult_dict = json.load(openfile)
        chain=item[-5]
        for index,row in new_df.iterrows():
            #my_mut=new_df.Mutation[index]
            #new_df.loc[index,'Mutation'] = my_mut[0]+chain+my_mut[1:]
            if chain in mult_dict.keys():
                new_df.loc[index, "ddG (kcal/mol)"] = mult_dict[chain]*new_df['ddG (kcal/mol)'][index]
                new_df.name=pdb_id
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
        DDG_col = df['ddG (kcal/mol)']
        ##use name_dict to get the correct name for the pdb file
        name_of_df=name_dict.get(df.name)
        combined_df[name_of_df] = DDG_col

    #making site column
    site=[]
    for index,row in combined_df.iterrows():
        site.append(row[0][2:-1])
    combined_df.insert(loc=0, column='site',value=site)

    #making unfolded column
    combined_df.insert(loc=5, column='unfolded',value=0)
    
    combined_df.to_csv(folder+'/'+prot_name+'_thermo_ddgs.csv', index=False)
    
    return combined_df