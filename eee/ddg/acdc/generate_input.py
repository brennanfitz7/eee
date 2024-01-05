"""
Functions to generate input for acdc_nn
"""

from eee.ddg.make_mutation_file import make_mutation_file
from eee.io import read_structure
from eee.io import write_fasta
from eee.core.data import AA_3TO1

import pandas as pd

import subprocess



def _make_psi_file(pdb_file:str, psi_file:str,hhblits_path:str,uniref_path:str,remove_multiple_models:bool):
    """
    Creates a psi file.

    Parameters
    ----------
    
    pdb_file : str
        pdb_file/path to pdb_file that will be used to create the fasta file
    
    psi_file : str
        desired name for the created psi file

    hhblits_path : str
        path to hhblits
    
    uniref_path : str
        path to uniref database

    remove_multiple_models : bool
        bool that determines whether multiple models are kept in the pdb file after read_structure

    Returns
    -------
    None
    """
    fasta_file=pdb_file.split('.')[0]+'.fasta' ##consider if later you want files to be 1abc.pdb_clean.fasta so you know they're from the clean version--something to consider

    pdb_df=read_structure(pdb_file,
                          remove_multiple_models=remove_multiple_models)

    
    write_fasta(pdb_df,fasta_file)
    
    verbose= True


    cmd=[hhblits_path,'-d',uniref_path,'-i',fasta_file,'-cpu','6','-n','2','-opsi',psi_file]


    popen = subprocess.Popen(cmd,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT,
                         universal_newlines=True)

    for line in popen.stdout:
        if verbose:
            print(line,end="",flush=True)
        
    # Check for success
    return_code = popen.wait()
    if return_code != 0:
        err = "Program failed.\n"
        raise RuntimeError(err)
        
        

def _make_prof_file_from_psi(psi_file:str, prof_file:str):
    """
    Creates a prof file from a psi file.

    Parameters
    ----------
    
    psi_file : str
        psi file for the protein
    
    prof_file : str
        prof file for the protein

    Returns
    -------
    None
    """
    verbose= True

    cmd=['ddgun','mkprof', psi_file]


    popen = subprocess.Popen(cmd,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT,
                         universal_newlines=True)
    
    f = open(prof_file, 'w')
    
    for line in popen.stdout:
        f.write(line)
    f.close()

        
    # Check for success
    return_code = popen.wait()
    if return_code != 0:
        err = "Program failed.\n"
        raise RuntimeError(err)

def _make_prof_file(pdb_file:str,psi_file:str,prof_file:str,hhblits_path:str,uniref_path:str,remove_multiple_models:bool):
    """
    Creates a prof file from a pdb file.

    Parameters
    ----------
    pdb_file : str
        pdb file name/path to pdb file

    psi_file : str
        psi file to be created
    
    prof_file : str
        prof file for the protein

    hhblits_path : str
        path to hhblits
    
    uniref_path : str
        path to uniref database

    remove_multiple_models : bool
        bool that determines whether multiple models are kept in the pdb file after read_structure

    Returns
    -------
    None
    """
    
    _make_psi_file(pdb_file=pdb_file, 
                   psi_file=psi_file, 
                   hhblits_path=hhblits_path, 
                   uniref_path=uniref_path,
                   remove_multiple_models=remove_multiple_models)

   
    _make_prof_file_from_psi(psi_file=psi_file, prof_file=prof_file)


def _acdc_nn_format(pdb_file:str, prof_file:str, tsv_file:str,remove_multiple_models:bool,just_a_test=True):
    """
    Makes a tab separated file in the format for ACDC_NN batch input.

    Parameters
    ----------
        
    pdb_file : str
        pdb file of protein of interest
        
    prof_file : str
        prof file to be created
    
    tsv_file : str
        tsv file to be created

    just_a_test : bool
        bool that determines whether the full mutation file is produced or only a test file of 10 mutations

    Returns
    -------
    None
    """

    ##mut_file=_make_mutation_file(pdb_file) change this when you're not running tests anymore


    mut_file=make_mutation_file(pdb_file,
                                remove_multiple_models=remove_multiple_models)


    #if just_a_test is True then only first 10 mutaitons will be in file
    if just_a_test==True:
        mut_file=mut_file.iloc[0:10]

    else:
        print("just_a_test argument entered incorrectly in acdc_nn_format")
    mutation_col=[]
    for index, row in mut_file.iterrows():
        mutation_col.append(AA_3TO1.get(row[1])+row[2]+AA_3TO1.get(row[3]))
    mut_file.drop(columns=['wt_res','res_pos','mut_res'], inplace=True)
    mut_file.insert(loc=0, column= 'SUB', value= mutation_col)
    mut_file.insert(loc=1, column= 'PROFILE', value=prof_file)
    mut_file.insert(loc=2, column= 'PDB', value=pdb_file)
    mut_file.rename(columns={"chain":"CHAIN"},inplace=True)
    
    mut_file.to_csv(tsv_file, sep="\t",header=False, index=False)
    
def generate_input(pdb_file:str, hhblits_path:str, uniref_path:str, remove_multiple_models:bool,just_a_test=True,):
    """
    Creates a prof file and an tab separated file in the format for ACDC_NN batch input. 

    Parameters
    ----------
        
    pdb_file : str
        pdb file of protein of interest

    fasta_file : str
        fasta file name/path to fasta file
    
    hhblits_path : str
        path to hhblits
    
    uniref_path : str
        path to uniref database 

    just_a_test : bool
        bool that determines whether the full mutation file is produced or only a test file of 10 mutations

    remove_multiple_models : bool
        bool that determines whether multiple models are kept in the pdb file after read_structure

    Returns
    -------
    None
    """ 
    
    pdb_id=pdb_file.split('.')[0]
    psi_file=pdb_id+'.psi'
    prof_file=pdb_id+'.prof'
    tsv_file=pdb_id+'_acdc_muts.tsv'


    _make_prof_file(pdb_file=pdb_file, 
                    psi_file=psi_file,
                    prof_file=prof_file, 
                    hhblits_path=hhblits_path, 
                    uniref_path=uniref_path, 
                    remove_multiple_models=remove_multiple_models) 


    _acdc_nn_format(pdb_file=pdb_file, 
                    prof_file=prof_file, 
                    tsv_file=tsv_file, 
                    just_a_test=just_a_test, 
                    remove_multiple_models=remove_multiple_models)
