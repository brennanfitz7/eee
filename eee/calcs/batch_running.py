from eee.analysis.dms_epistasis import dms_epistasis

import pandas as pd
import subprocess

def batch_run_dms_epi(on_off_df,prot_name:str):
    """
    must be run in a file with an ensemble.csv, conditions.csv, simulation.json, ddg.csv.
    """
    
    conditions_df=pd.read_csv('conditions.csv')
    for idx,row in on_off_df.iterrows():
        conditions_df.at[0, 'oht'] = on_off_df['off'][idx]
        conditions_df.at[1,'oht']=on_off_df['on'][idx]
        conditions_df.to_csv('conditions.csv')

        cmd=['eee-run-calculation','simulation.json',prot_name+"_on_"+on_off_df['on'][idx]+'_off_'+on_off_df['off'][idx]]

        verbose=True

        #using call instead of popen so that the whole thing is run before conditions.csv is rewritten
        call= subprocess.call(cmd,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             universal_newlines=True)

        for line in call.stdout:
            if verbose:
                print(line,end="",flush=True)

        # Check for success
        return_code = call.wait()
        if return_code != 0:
            err = "Program failed.\n"
            raise RuntimeError(err)


        epi_df=dms_epistasis(prot_name+"_on_"+on_off_df['on'][idx]+'_off_'+on_off_df['off'][idx])
        epi_df.to_csv(prot_name+"_on_"+on_off_df['on'][idx]+'_off_'+on_off_df['off'][idx]+'_epi.csv')
