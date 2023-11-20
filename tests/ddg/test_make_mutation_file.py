from eee.io.read_structure import read_structure
from eee.ddg.make_mutation_file import make_mutation_file

import pandas as pd


def test__make_mutation_file(test_pdbs):
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
        result_df=make_mutation_file(file)
        assert issubclass(type(result_df),pd.DataFrame)