#from eee.analysis.dms_epistasis import dms_epistasis

import pandas as pd
import shutil
import os

from eee.calcs import read_json

def batch_run_dms_epi(json_file,
                      midpoints:list,
                      switch_range:list,
                      ligand:str,
                      prot_name:str,
                      on_is_higher:bool,
                      use_stored_seed=False,
                      overwrite=False):
    """
    must be run in a file with an ensemble.csv, conditions.csv, simulation.json, ddg.csv.

    on_is_higher indicates that the observable should be "on" at the highest point in the range. If on_is_higher is False, the observable should be "off"
    """
    on=[]
    off=[]
    for point in midpoints:
        for number in switch_range:
            if on_is_higher==True:
                on.append(point+number)
                off.append(point-number)
            if on_is_higher==False:
                on.append(point-number)
                off.append(point+number)

    on_off_df=pd.DataFrame(data=({'on':on,'off':off}))

    conditions_df=pd.read_csv('conditions.csv')
    for idx,row in on_off_df.iterrows():
        on=on_off_df['on'][idx]
        off=on_off_df['off'][idx]
        conditions_df.loc[conditions_df.fitness_fcn=='on', ligand] = on
        conditions_df.loc[conditions_df.fitness_fcn=='on', ligand] = on
        conditions_df.to_csv('conditions.csv',index=False)

        out_dir=prot_name+"_on_"+str(on)+'_off_'+str(off)

        es, kwargs = read_json(json_file=json_file,
                        use_stored_seed=use_stored_seed)
        
        if os.path.exists(out_dir):
            if overwrite:
                shutil.rmtree(out_dir)
            else:
                err = f"\noutput_directory '{out_dir}' already exists\n\n"
                raise FileExistsError(err)

        print()
        print(es.get_calc_description(kwargs))
        print(f"\nWriting result to '{out_dir}' directory\n",flush=True)
            
        es.run(output_directory=out_dir,**kwargs)

        #cannot run dms_epistasis within eee because it meses with importing
        #epi_df=dms_epistasis(prot_name+"_on_"+on_off_df['on'][idx]+'_off_'+on_off_df['off'][idx])
        #epi_df.to_csv(prot_name+"_on_"+on_off_df['on'][idx]+'_off_'+on_off_df['off'][idx]+'_epi.csv')
