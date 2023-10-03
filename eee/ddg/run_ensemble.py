

from eee.ddg import acdc
from eee.ddg import foldx
from eee.structure import sync_structures

import pandas as pd

import glob
import sys
import argparse
import shutil



def run_ensemble(pdb_csv:str,prot_name:str,module:str,just_a_test=True):
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
    
    just_a_test : bool
        bool that determines whether the full mutation file is run or only a test file of 10 mutations
        
    Returns
    -------
    None
    """
    pdb_df=pd.read_csv(pdb_csv)
    pdb_list=pdb_df['PDB']
    name_list=pdb_df['NAME']
    name_dict = dict(zip(pdb_list, name_list)) #sets pdbs as keys and names as values
    module_dict={'acdc':acdc,'foldx':foldx}
    calculator=module_dict.get(module)
    if calculator==acdc:
        acdc_info=pdb_df['ACDC'].tolist()
        hhblits_path=acdc_info[0]
        uniref_path=acdc_info[1]
    
    sync_structures(structure_files=pdb_list, out_dir=prot_name)

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
            if just_a_test==True:
                calculator.generate_input(pdb_file)
            elif just_a_test==False:
                calculator.generate_input(pdb_file, just_a_test=False)
            else:
                print('just_a_test argument entered incorrectly in run_ensemble')
            calculator.ddg_calc(muts_file=pdb_id+'_'+calculator+'_muts.tsv', pdb_file=pdb_file)
            calculator.convert_to_df('PS'+pdb_file[0:-4]+'_scanning_output.txt')
            

            
        if calculator==acdc:
            pdb_file=str(prot_name+'/'+pdb)
            pdb_id=pdb_file.split('.')[0]
            if just_a_test==True:
                calculator.generate_input(pdb_file,hhblits_path=hhblits_path,uniref_path=uniref_path)
            elif just_a_test==False:
                calculator.generate_input(pdb_file,hhblits_path=hhblits_path, uniref_path=uniref_path, just_a_test=False)
            else:
                print('just_a_test argument entered incorrectly in run_ensemble')
            #run ddg_calc and convert_to_df
            calculator.ddg_calc(pdb_id+'_'+calculator+'_muts.tsv')
            ddg_df=calculator.convert_to_df(pdb_id+'_'+calculator+'_raw_ddgs.txt')


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

    mut_source=clean_df_list[0]
    combined_df=mut_source[['Mutation']].copy()

    for df in clean_df_list:
        DDG_col = df['DDG']
        ##use name_dict to get the correct name for the pdb file
        name_of_df=name_dict.get(df.name)
        combined_df[name_of_df] = DDG_col
    
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
    
