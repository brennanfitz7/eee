from eee.ddg.acdc.convert_to_df import convert_to_df ##cannot test fully without example documents--ask Mike
from eee.data import AA_3TO1
from eee.io import read_structure


import pandas as pd
import os
from eee.io import read_structure

#need example outputs from acdc to test this function
def test_convert_to_df(test_acdc_ddg_outputs):
    for ddg_output in test_acdc_ddg_outputs:
        pdb_id=ddg_output.split('_')[0]
        muts_file=pdb_id+'_acdc_muts.tsv'
        output_df=pdb_id+'_acdc_ddg_df.csv'

        with open(ddg_output,'r') as ddgs:
            output_list = ddgs.read().splitlines()
        
        ddg_list=[]
        for item in output_list:
            if item[3:len(item)].isdigit()==True:
                ddg_list.append(item)
            else:
                continue
        
        #test that everything that ended up in ddg list is a number
        for item in ddg_list:
            num=float(item)
            assert issubclass(type(num),float)

        input_df=pd.read_table(muts_file, delimiter='\t',names=['Mutation', 'Prof', 'PDB','Chain'])
        
        ddg_df = pd.DataFrame(
        {'Mutation': input_df['Mutation'],'DDG': ddg_list})

        #test that these the columns exist in the dataframe and that they have the correct lengths
        assert issubclass(type(ddg_df),pd.DataFrame)
        assert len(ddg_df.columns)==2
        assert len(ddg_df['Mutation'])==len(input_df['Mutation'])
        assert len(ddg_df['DDG'])==len(ddg_list)
        assert len(ddg_df['Mutation'])==len(ddg_df['DDG'])

        ddg_df.to_csv(output_df, index=False)

        #check ddg_df was produced to csv
        assert os.path.isfile(output_df)==True
    