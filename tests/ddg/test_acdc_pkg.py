from eee.ddg.acdc_pkg import _make_psi_file ##cannot test on personal computer
from eee.ddg.acdc_pkg import _make_prof_file_from_psi ##cannot test on personal computer
from eee.ddg.acdc_pkg import _make_prof_file ##cannot test on personal computer
from eee.ddg.acdc_pkg import _make_mutation_file
from eee.ddg.acdc_pkg import _acdc_nn_format
from eee.ddg.acdc_pkg import generate_input ##cannot test on personal computer
from eee.ddg.acdc_pkg import ddg_calc ##cannot test on personal computer
from eee.ddg.acdc_pkg import convert_to_df ##cannot test fully without example documents--ask Mike



import sys
import argparse
import pandas as pd
from eee.io import read_structure

AA = [("A","ALA"),
      ("C","CYS"),
      ("D","ASP"),
      ("E","GLU"),
      ("F","PHE"),
      ("G","GLY"),
      ("H","HIS"),
      ("I","ILE"),
      ("K","LYS"),
      ("L","LEU"),
      ("M","MET"),
      ("N","ASN"),
      ("P","PRO"),
      ("Q","GLN"),
      ("R","ARG"),
      ("S","SER"),
      ("T","THR"),
      ("V","VAL"),
      ("W","TRP"),
      ("Y","TYR")]

AA_3TO1 = dict([(a[1],a[0]) for a in AA])
AA_1TO3 = dict([(a[0],a[1]) for a in AA])

def test__make_mutation_file(test_pdbs:str):
    for file in test_pdbs:
        raw_prot=read_structure(file)
        assert issubclass(type(raw_prot),pd.DataFrame)

        prot_seq=raw_prot[raw_prot['class'] == "ATOM"].loc[:,["chain",'resid','resid_num']].drop_duplicates(subset=["resid_num","chain"], keep='first')
        assert len(prot_seq.columns)==3
        assert len(prot_seq.index)<(len(raw_prot.index)/4)
    
        mutation_dict={'chain':[],'wt_res':[],'res_pos':[],'mut_res':[]}
        assert issubclass(type(mutation_dict),dict)
         
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
        #test that each residue has 19 possible mutations in this dictionary
        assert len(mutation_dict["chain"])==(len(prot_seq.index)*19)
        assert len(mutation_dict["wt_res"])==(len(prot_seq.index)*19)
        assert len(mutation_dict["res_pos"])==(len(prot_seq.index)*19)
        assert len(mutation_dict["mut_res"])==(len(prot_seq.index)*19)
    
        #test that the mutation dataframe created is indeed a dataframe
        mutation_df=pd.DataFrame.from_dict(mutation_dict)
        assert issubclass(type(mutation_df),pd.DataFrame)

        #test that ultimately _make_mutation_df does return a dataframe
        result_df=_make_mutation_file(file)
        assert issubclass(type(result_df),pd.DataFrame)


def test__acdc_nn_format(test_pdbs):
    for file in test_pdbs:
        mut_file=_make_mutation_file(file).iloc[0:10]
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

    
        #mut_file.to_csv('tsv_file', sep="\t",header=False, index=False)
#this test is not done!
def test_convert_to_df(test_pdbs):
    for p in test_pdbs:
        assert issubclass(type(convert_to_df(p)),pd.DataFrame)