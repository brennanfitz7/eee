from eee.ddg.foldx.ddg_calc import ddg_calc

import pandas as pd
import os

def test_ddg_calc(test_pdbs,test_muts):
    """
    Note to me: need to make csv with test_pdbs and test mut files
    (each row has one pdb and one mut_file)
    """
    
    for pdb in test_pdbs:
        for mut_file in test_muts:
            if mut_file[0:4]==pdb[0:4]:
                ddg_calc(muts_file=mut_file, pdb_file=pdb)
                
                #test that file was produced
                assert os.path.isfile('PS_'+pdb[0:-4]+'_scanning_output.txt')