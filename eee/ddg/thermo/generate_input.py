from eee.structure.sync_structures import sync_structures
from eee.structure.chains_and_mults import chains_and_multipliers

import glob

def generate_input(list_of_pdbs:list, prot_name:str, all_models_necessary:bool):
    
    new_prot_name=prot_name+'_synced'
    sync_structures(list_of_pdbs, 
                    out_dir=new_prot_name, 
                    all_models_necessary=all_models_necessary)
    
    synced_pdbs=glob.glob(prot_name+'/*.pdb')

    for pdb in synced_pdbs:
        chains_and_multipliers(pdb,save_dict=True)