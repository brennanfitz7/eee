#from eee.analysis.dms_epistasis import dms_epistasis

import pandas as pd
import subprocess

def batch_run_dms_epi(on_off_df,prot_name:str):
    """
    must be run in a file with an ensemble.csv, conditions.csv, simulation.json, ddg.csv.
    """
    
    conditions_df=pd.read_csv('conditions.csv')
    for idx,row in on_off_df.iterrows():
        on=on_off_df['on'][idx]
        off=on_off_df['off'][idx]
        conditions_df.at[0, 'oht'] = off
        conditions_df.at[1,'oht']=on
        conditions_df.to_csv('conditions.csv')

        out_dir=prot_name+"_on_"+on+'_off_'+off

        cmd=['eee-run-calculation','simulation.json',out_dir]

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

        #cannot run dms_epistasis within eee because it meses with importing
        #epi_df=dms_epistasis(prot_name+"_on_"+on_off_df['on'][idx]+'_off_'+on_off_df['off'][idx])
        #epi_df.to_csv(prot_name+"_on_"+on_off_df['on'][idx]+'_off_'+on_off_df['off'][idx]+'_epi.csv')
