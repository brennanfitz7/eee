"""
Functions to generate input for acdc_nn
"""


from eee.io import read_structure
from eee.io import write_fasta
from eee.data import AA_3TO1

import pandas as pd

import subprocess



def _make_psi_file(pdb_file:str, psi_file:str):
    """
    Creates a psi file.

    Parameters
    ----------
    
    pdb_file : str
        pdb_file/path to pdb_file that will be used to create the fasta file
    
    psi_file : str
        desired name for the created psi file

    Returns
    -------
    None
    """
    fasta_file=pdb_file.split('.')[0] ##consider if later you want files to be 1abc.pdb_clean.fasta so you know they're from the clean version--something to consider
    pdb_df=read_structure(pdb_file)
    write_fasta(pdb_df,fasta_file)
    
    verbose= True


    #both of these paths have to be changed to make it work for anyone outside of my specfic spock account
    cmd=['/home/brennanfitz7/miniconda3/bin/hhblits','-d','/home/brennanfitz7/ACDC_NN/UniRef30_2023_02/UniRef30_2023_02','-i',fasta_file,'-cpu','6','-n','2','-opsi',psi_file]


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

def _make_prof_file(pdb_file:str,psi_file:str,prof_file:str):
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

    Returns
    -------
    None
    """
    _make_psi_file(pdb_file=pdb_file, psi_file=psi_file )
    _make_prof_file_from_psi(psi_file=psi_file, prof_file=prof_file)


def _make_mutation_file(pdb_file:str):
    """
    Makes a dataframe with all possible mutations from a pdb file. 

    Parameters
    ----------
        
    pdb_file : str
        pdb file of protein of interest

    Returns
    -------
    Dataframe with all possible mutations from a pdb file. 
    """
    
    #get dataframe from pdb, select for atoms, and drop duplicates
    raw_prot=read_structure(pdb_file)
    prot_seq=raw_prot[raw_prot['class'] == "ATOM"].loc[:,["chain",'resid','resid_num']].drop_duplicates(subset=["resid_num","chain"], keep='first')
    
    #create new mutation dictionary
    mutation_dict={'chain':[],'wt_res':[],'res_pos':[],'mut_res':[]}
    
    aa_list=["ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR"]
    
    #populate mutation dictionary with info from original df
    for index, row in prot_seq.iterrows():
        for aa in aa_list:
            if row[1]!=aa:
                mutation_dict["chain"].append(row[0])
                mutation_dict["wt_res"].append(row[1])
                mutation_dict["res_pos"].append(row[2])
                mutation_dict["mut_res"].append(aa)
            else:
                continue
    
    #create a dataframe from the mutation dictionary
    mutation_df=pd.DataFrame.from_dict(mutation_dict)
    
    return mutation_df


def _acdc_nn_format(pdb_file:str, prof_file:str, tsv_file:str):
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

    Returns
    -------
    None
    """

    ##mut_file=_make_mutation_file(pdb_file) change this when you're not running tests anymore
    mut_file=_make_mutation_file(pdb_file).iloc[0:10]
    mutation_col=[]
    for index, row in mut_file.iterrows():
        mutation_col.append(AA_3TO1.get(row[1])+row[2]+AA_3TO1.get(row[3]))
    mut_file.drop(columns=['wt_res','res_pos','mut_res'], inplace=True)
    mut_file.insert(loc=0, column= 'SUB', value= mutation_col)
    mut_file.insert(loc=1, column= 'PROFILE', value=prof_file)
    mut_file.insert(loc=2, column= 'PDB', value=pdb_file)
    mut_file.rename(columns={"chain":"CHAIN"},inplace=True)
    
    mut_file.to_csv(tsv_file, sep="\t",header=False, index=False)
    
def generate_input(pdb_file:str):
    """
    Creates a prof file and an tab separated file in the format for ACDC_NN batch input. 

    Parameters
    ----------
        
    pdb_file : str
        pdb file of protein of interest

    fasta_file : str
        fasta file name/path to fasta file

    Returns
    -------
    None
    """ 
    
    pdb_id=pdb_file.split('.')[0]
    psi_file=pdb_id+'.psi'
    prof_file=pdb_id+'.prof'
    tsv_file=pdb_id+'_ddg_input.tsv'
    

    _make_prof_file(pdb_file=pdb_file, psi_file=psi_file,prof_file=prof_file) 
    _acdc_nn_format(pdb_file=pdb_file, prof_file=prof_file, tsv_file=tsv_file)