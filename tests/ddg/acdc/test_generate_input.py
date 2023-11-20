from eee.ddg.acdc.generate_input import _make_psi_file 
from eee.ddg.acdc.generate_input import _make_prof_file_from_psi 
from eee.ddg.acdc.generate_input import _make_prof_file 
from eee.ddg.make_mutation_file import make_mutation_file
from eee.ddg.acdc.generate_input import _acdc_nn_format
from eee.data import AA_3TO1


import pandas as pd
import os



def test__make_psi_file(test_acdc_csv):

    """
    Note to me: test_acdc_csv should contain columns labeled PDB, HHBLITS, and UNIREF
    PDB will have a list of all test_pdb's. HHBLITS and UNIREF will only have the path.
    """
    csv_df=pd.read_csv(test_acdc_csv)
    pdb_list=csv_df['PDB'].tolist()
    hhblits_path=csv_df['HHBLITS'][0]
    uniref_path=csv_df['UNIREF'][0]
    for pdb in pdb_list:
        pdb_id=pdb.split('.')[0]
        psi_file=pdb_id+'.psi'
        tsv_file=pdb_id+'_acdc_muts.tsv'

        _make_psi_file(pdb_file=pdb, psi_file=psi_file, hhblits_path=hhblits_path, uniref_path=uniref_path)

        #test that each psi file was made
        assert os.path.isfile(psi_file)==True

def test__make_prof_file_from_psi(test_psi_files):
    for psi_file in test_psi_files:
        prof_file=psi_file.split('.')[0]+'.prof'
        _make_prof_file_from_psi(psi_file, prof_file)

        #test that each prof file exists
        assert os.path.isfile(prof_file)==True

def test__make_prof_file(test_acdc_csv):
    """
    Note to me: test_acdc_csv should contain columns labeled PDB, HHBLITS, and UNIREF
    PDB will have a list of all test_pdb's. HHBLITS and UNIREF will only have the path.
    """
    csv_df=pd.read_csv(test_acdc_csv)
    pdb_list=csv_df['PDB'].tolist()
    hhblits_path=csv_df['HHBLITS'][0]
    uniref_path=csv_df['UNIREF'][0] 

    for pdb in pdb_list:
        pdb_id=pdb.split('.')[0]
        psi_file=pdb_id+'.psi'
        prof_file=pdb_id+'.prof'
        _make_prof_file(psi_file=psi_file, prof_file=prof_file, hhblits_path=hhblits_path, uniref_path=uniref_path)

        #test that each psi file was made
        assert os.path.isfile(psi_file)==True

        #test that each prof file exists
        assert os.path.isfile(prof_file)==True


def test__acdc_nn_format(test_pdbs):
    for file in test_pdbs:
        mut_file=make_mutation_file(file).iloc[0:10]
        mutation_col=[]
        for index, row in mut_file.iterrows():
            mutation_col.append(AA_3TO1.get(row[1])+row[2]+AA_3TO1.get(row[3]))
        #make sure mutation_col was populated correctly
        assert len(mut_file.index)==len(mutation_col)
        mut_file.drop(columns=['wt_res','res_pos','mut_res'], inplace=True)
        #make sure there is only one column
        assert len(mut_file.columns)==1
        mut_file.insert(loc=0, column= 'SUB', value= mutation_col)
        mut_file.insert(loc=1, column= 'PROFILE', value='prof_file')
        mut_file.insert(loc=2, column= 'PDB', value=file)
        mut_file.rename(columns={"chain":"CHAIN"},inplace=True)
        #assert that there are four columns
        assert len(mut_file.columns)==4
        #assert that "chain" is no longer a column
        assert "chain" not in mut_file
        
        pdb_id=file.split('.')[0]
        tsv_file=pdb_id+'_acdc_muts.tsv'

        mut_file.to_csv(tsv_file, sep="\t",header=False, index=False)

        #assert that tsv file was made
        assert os.path.isfile(tsv_file)==True

def test_generate_input(test_acdc_csv):
    """
    Note to me: test_acdc_csv should contain columns labeled PDB, HHBLITS, and UNIREF
    PDB will have a list of all test_pdb's. HHBLITS and UNIREF will only have the path.
    """
    csv_df=pd.read_csv(test_acdc_csv)
    pdb_list=csv_df['PDB'].tolist()
    hhblits_path=csv_df['HHBLITS'][0]
    uniref_path=csv_df['UNIREF'][0] 

    for pdb in pdb_list:
        pdb_id=pdb.split('.')[0]
        psi_file=pdb_id+'.psi'
        prof_file=pdb_id+'.prof'
        tsv_file=pdb_id+'_acdc_muts.tsv'

        _make_prof_file(pdb_file=pdb,psi_file=psi_file,prof_file=prof_file)
        #assert that prof_file was made (if prof_file was made,can assume psi_file was as well)
        assert os.path.isfile(prof_file)==True

        _acdc_nn_format(pdb_file=pdb, prof_file=prof_file, tsv_file=tsv_file)
        assert os.path.isfile(tsv_file)==True
