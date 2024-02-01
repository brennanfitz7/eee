

from eee.ddg import acdc
from eee.ddg import foldx
from eee.structure import sync_structures

import pandas as pd

import glob
import sys
import argparse
import shutil
import json


def regularize_data(df, stats_dict,stdev_rosetta=13.228994023873321):
    """
    Regularize data from a different ∆∆G calculator to Rosetta's standard

    Parameters
    ----------
    df : pandas DataFrame
        dataframe output from convert_to_df function

    stats_dict : dict
        dictionary with mean and stdev for whichever calculator the data is from. Should contain keys "mean" and "stdev"
    
    stdev_rosetta : float
        the standard deviation for rosseta data--will be updated as more data is obtained
        
    Returns
    -------
    DataFrame with regularized data
    """
    mean_other=stats_dict.get('mean')
    stdev_other=stats_dict.get('stdev')
    for column in df.columns[2:-1]:
        df[column] = df[column].map(lambda item: ((((item-mean_other)/stdev_other)*stdev_rosetta)+((stdev_rosetta/stdev_other)*mean_other)))

    return df

def run_ensemble(pdb_csv:str,prot_name:str,module:str,all_models_necessary:bool,just_a_test=True):
    """
    Runs three pdb file (cif files eventually) and create ddg csv.

    Parameters
    ----------
    pdb_csv : str
        file name or path to file of csv file. File should contain columns NAME (column contains name designations, such as 'apo') and PDB (pdb file or filepath)
        if calculator is acdc, should also contain column ACDC with hhblits path and uniref path, in that order

    prot_name : str
        the name of the protein. This will determine the name of the output file created by sync structures.
    
    module : str
        the name of the module being used to calculate DDG. For now, must be 'acdc' or 'foldx'.

    all_models_necessary: bool, default=True
        bool for whether all models in a pdb are necessary for analysis. Multiple models should be removed for NMR structures but kept for biological assemblies that contain multiple asymmetric units. 
    
    just_a_test : bool
        bool that determines whether the full mutation file is run or only a test file of 10 mutations
        
    Returns
    -------
    Dataframe with ∆∆G's for each mutation for each residue for every file in the ensemble
    """
    
    pdb_df=pd.read_csv(pdb_csv)
    pdb_list=pdb_df['PDB']
    name_list=pdb_df['NAME']
    name_dict = dict(zip(pdb_list, name_list)) #sets pdbs as keys and names as values
    #module dict to take it from string to actual module
    module_dict={'acdc':acdc,'foldx':foldx}
    calculator=module_dict.get(module)
    #dicts for regularizing data
    foldx_dict={'mean': 4.043667500801383, 'stdev': 28.70340193260126}
    acdc_dict={'mean': 0.24449510098968205, 'stdev': 0.6223300393521769}
    if calculator==acdc:
        acdc_info=pdb_df['ACDC'].tolist()
        hhblits_path=acdc_info[0]
        uniref_path=acdc_info[1]
    
    sync_structures(structure_files=pdb_list, out_dir=prot_name, all_models_necessary=all_models_necessary)
    #this is the point at which models have either been incorporated or thrown out

    synced_pdbs=glob.glob('*.pdb',root_dir=prot_name)
    #this returns a list of pdb files without the file path

    
    #for foldx, need to move files into whatever file I'm running things in
    if calculator==foldx:
        for pdb in synced_pdbs: 
            shutil.move(prot_name+'/'+pdb,pdb)


    df_list=[]

    for pdb in synced_pdbs:
        if calculator==foldx:
            pdb_file=pdb
            pdb_id=pdb_file.split('.')[0]
            calculator.generate_input(pdb_file,
                                      just_a_test=just_a_test)
        
                
            calculator.ddg_calc(muts_file=pdb_id+'_'+module+'_muts.txt', pdb_file=pdb_file)
            ddg_df=calculator.convert_to_df('PS_'+pdb_file[0:-4]+'_scanning_output.txt')
            #This will regularize the data to rosetta's standard
            #ddg_df=regularize_data(ddg_df, foldx_dict)
            #^^MAY ADD THIS BACK LATER
            #multiply each DDG by the multiplier provided by the pickle (based on the number of identical chains)
            with open(pdb_id+'_ddg_mult.json', 'r') as openfile:
                mult_dict = json.load(openfile)
            for index,row in ddg_df.iterrows():
                chain=ddg_df['Mutation'][index][1]
                if chain in mult_dict.keys():
                    ddg_df.loc[index, "DDG"] = mult_dict[chain]*ddg_df['DDG'][index]
            ddg_df.name=pdb.split('_')[0]
            df_list.append(ddg_df)

            
        if calculator==acdc:
            pdb_file=str(prot_name+'/'+pdb)
            pdb_id=pdb_file.split('.')[0]
            calculator.generate_input(pdb_file,
                                      hhblits_path=hhblits_path,
                                      uniref_path=uniref_path,
                                      just_a_test=just_a_test)
            #run ddg_calc and convert_to_df
            calculator.ddg_calc(pdb_id+'_'+module+'_muts.tsv')
            ddg_df=calculator.convert_to_df(pdb_id+'_'+module+'_raw_ddgs.txt')
            #This will regularize the data to rosetta's standard
            #ddg_df=regularize_data(ddg_df, acdc_dict)
            #^^MAY ADD THIS BACK LATER
            #multiply each DDG by the multiplier provided by the pickle (based on the number of identical chains)
            with open(pdb_id+'_ddg_mult.json', 'r') as openfile:
                mult_dict = json.load(openfile)
            for index,row in ddg_df.iterrows():
                chain=ddg_df['Mutation'][index][1]
                if chain in mult_dict.keys():
                    ddg_df.loc[index, "DDG"] = mult_dict[chain]*ddg_df['DDG'][index]
            ddg_df.name=pdb.split('_')[0]
            df_list.append(ddg_df)


    #create list of sets of each df's mutation column 
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
        DDG_col = df['DDG']
        ##use name_dict to get the correct name for the pdb file
        name_of_df=name_dict.get(df.name)
        combined_df[name_of_df] = DDG_col

    #making site column
    site=[]
    for index,row in combined_df.iterrows():
        site.append(row[0][1:-1])
    combined_df.insert(loc=0, column='site',value=site)

    #making unfolded column
    combined_df.insert(loc=5, column='unfolded',value=0)
    
    combined_df.to_csv(prot_name+'_'+module+'_ddgs.csv', index=False)
    
    return combined_df
    
    
    
def main(argv=None):
    
    if argv is None:
        argv = sys.argv[1:]
        

    parser = argparse.ArgumentParser(description='Find the DDG for each mutation in an ensemble')
    parser.add_argument('-p','--pdb_csv',metavar='', required=True, help='CSV file with pdb names and files')
    parser.add_argument('-n','--prot_name',metavar='', required=True, help='Name of protein')
    parser.add_argument('-m','--module',metavar='', required=True, help='Module to be used for DDG calculations')
    args = parser.parse_args(argv)

    
    run_ensemble(pdb_csv=args.pdb_csv, 
                 prot_name=args.prot_name,
                 module=args.module)
        

if __name__ == "__main__":
    main()
    
