
from eee.ddg.make_mutation_file import make_mutation_file
from eee.data import AA_3TO1

import os

def test_generate_input(test_pdbs):
    for pdb in test_pdbs:
        mut_file_long=make_mutation_file(pdb)

        #making a list with all the different mutations    
        mut_file=mut_file_long.drop_duplicates(subset=["res_pos"], keep='first')
        assert len(mut_file_long.index)==len(mut_file.index)*19

        mutation_col=[]
        for index, row in mut_file.iterrows():
            mutation_col.append(AA_3TO1.get(row[1])+row[0]+row[2]+'a')

        #test that all the mutations from mut made it into mutation_col
        assert len(mutation_col)==len(mut_file.index)

        mut_string=','.join(mut for mut in mutation_col)

        #test that the length makes sense (assuming 4 characters per mutation +1 comma per each and no comma on end)
        assert len(mut_string)==(len(mutation_col)*5)-1

        txt_file=pdb.split('.')[0]+'_foldx_muts.txt'
        mut_txt = open(txt_file, "w")
        mut_txt.write(mut_string)
        mut_txt.close()

        #check that file was made
        os.path.isfile(txt_file)
