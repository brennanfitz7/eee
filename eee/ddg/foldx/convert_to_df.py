from eee.data import AA_3TO1

import pandas as pd

def convert_to_df(pdb_file:str):
#figure out how to add header without deleting first row--it's okay for now bc the first line is mutating the residue to itself
    #foldx is able to identify the protein tag alone--no path
    pdb_id=pdb_file.split('/')[-1].split('.')[0]
    #establishing the name of the output file to find--should be where I ran it regardless of where pdb files are
    ddg_txt='PS_'+pdb_id+'_scanning_output.txt'
    column_names=['Mutation','DDG']
    ddg_df=pd.read_csv(ddg_txt, sep='\t')

    ddg_df.set_axis(column_names, axis=1, inplace=True)
    for index,row in ddg_df.iterrows():
        one_letter_code=AA_3TO1.get(row[0][0:3])
        if one_letter_code!=row[0][5] and row[0][5].isupper():
            #replacing the old mutation format with a new one
            ddg_df.replace(row[0],one_letter_code+row[0][4]+row[0][5],inplace=True)
        else:
            ddg_df.drop(index,inplace=True)

    output_df=pdb_id+'_ddg_df.csv'
    ddg_df.to_csv(output_df, index=False)
    
    return ddg_df