

from eee.ddg import acdc
from eee.structure import sync_structures

import pandas as pd

import glob
import sys
import argparse
import shutil



def run_ensemble(pdb_csv:str,prot_name:str,module:str):
    """
    Runs three pdb file (cif files eventually) and create ddg csv.

    Parameters
    ----------
    pdb_csv : str
        file name or path to file of csv file. File should contain columns NAME (column contains name designations, such as 'apo') and PDB (pdb file or filepath)

    prot_name : str
        the name of the protein. This will determine the name of the output file created by sync structures.
    
    module : str
        the name of the module being used to calculate DDG. For now, must be 'acdc'.
        
    Returns
    -------
    None
    """
    pdb_df=pd.read_csv(pdb_csv)
    pdb_list=pdb_df['PDB']
    name_list=pdb_df['NAME']
    name_dict = dict(zip(pdb_list, name_list)) #sets pdbs as keys and names as values
    module_dict={'acdc':acdc}
    calculator=module_dict.get(module)
    
    sync_structures(structure_files=pdb_list, out_dir=prot_name)

    synced_pdbs=glob.glob('*.pdb',root_dir=prot_name)
    #this returns a list of pdb files without the file path
    #FIX THIS SO PDB's GO IN IN ORDER THEY ARE MADE NOT IN ALPHABETICAL ORDER
    #or maybe make a dictionary with the csv file--seems possible to create a dictionary from two lists
    
    df_list=[]
        
    for pdb in synced_pdbs:
        pdb_file=str(prot_name+'/'+pdb)
        calculator.generate_input(pdb_file)
        calculator.ddg_calc(pdb_file)
        pdb=calculator.convert_to_df(pdb_file)
        df.name=pdb.split('_')[0]
        df_list.append(pdb)

    
    ##make first column of dataframe with mutations from one of the files--they should all be the same
    ##maybe write code that checks that these are all the same? Is this necessary after sync_structures?
    #this function should not require running sync structures itself
    mut_source=df_list[0]
    combined_df=mut_source[['Mutation']].copy()

    for df in df_list:
        DDG_col = df['DDG']
        ##use name_dict to get the correct name for the pdb file
        name_of_df=name_dict.get(df.name)
        combined_df[name_of_df] = DDG_col
    
    combined_df.to_csv(prot_name+'/'+prot_name+'_combined_df.csv', index=False)
    
    
    
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
    
