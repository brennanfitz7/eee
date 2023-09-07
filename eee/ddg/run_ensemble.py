

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
    module_dict={'acdc':acdc}
    calculator=module_dict.get(module)
    
    sync_structures(structure_files=pdb_list, out_dir=prot_name)

    synced_pdbs=glob.glob('*.pdb',root_dir=prot_name)
    #this returns a list of pdb files without the file path
    
    df_list=[]
    
    if calculator!=acdc:
        for pdb in synced_pdbs:
            pdb_file=str(prot_name+'/'+pdb)
            calculator.generate_input(pdb_file)
            calculator.ddg_calc(pdb_file)
            pdb=calculator.convert_to_df(pdb_file)
            df_list.append(pdb)
           
            
    
    elif calculator==acdc:
        output=glob.glob('calculation_*/tmp-align_output.fasta')[0]
        with open(output,'r') as seqs:
            seq_list = seqs.read().split('>')
        c=0
        for seq in seq_list:
            simple_seq=seq.replace('\n','')[4:]
            #no placeholder hyphens
            ready_seq=simple_seq.replace('-','')
            pdb_tag=pdb_list[c].split('.')[0]
            #making fasta_seq
            fasta_seq='>'+pdb_tag+'\n'+ready_seq
            #write fasta_seq into fasta file
    
            file = open(pdb_tag+".fasta", "w")
            file.write(fasta_seq)
            file.close()
            c=c+1

        d=0 
        for pdb in synced_pdbs:  
            fasta_file=pdb_list[d].split('.')[0]+".fasta"
            pdb_file=str(prot_name+'/'+pdb)
            calculator.generate_input(pdb_file, fasta_file)
            calculator.ddg_calc(pdb_file)
            pdb=calculator.convert_to_df(pdb_file)
            df_list.append(pdb)
            d=d+1
            
    
    else:
        return "This should not happen."

    #this removes the calculation file containing the fasta alignments    
    shutil.rmtree(glob.glob('calculation_*')[0])
    
    ##make first column of dataframe with mutations from one of the files--they should all be the same
    ##maybe write code that checks that these are all the same? Is this necessary after sync_structures?
    mut_source=df_list[0]
    combined_df=mut_source[['Mutation']].copy()

    i=0
    for df in df_list:
        DDG_col = df['DDG']
        df.name=name_list[i]
        combined_df[df.name] = DDG_col
        i=i+1
    
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
    
