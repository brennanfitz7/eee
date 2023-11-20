from eee.ddg.acdc.ddg_calc import ddg_calc 

import os

def test_ddg_calc(test_acdc_muts_files):
    for muts_file in test_acdc_muts_files:
        ddg_calc(muts_file)
        #test that the proper file was produced
        os.path.isfile(muts_file.split('_')[0]+'_acdc_raw_ddgs.txt')

