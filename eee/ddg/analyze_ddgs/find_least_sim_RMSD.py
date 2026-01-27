from eee._private import logger

import pandas as pd
import numpy as np

import subprocess
import shutil
import glob
import os 


def find_longest_chain(chains_file):
    
    chains=[]
    
    with open(chains_file, 'r') as file:
        for line in file:
            chains.append(line.strip())
            
    if len(chains)>1:
        
        ens_dir=chains_file.split('_chains.txt')[0]
        
        chain_length_dict={}
        
        for chain in chains:
            
            length_chain=[]
        
            for ddg_file in glob.glob(ens_dir+'/*thermo_ddg_'+chain+'.csv'):
                
                df=pd.read_csv(ddg_file)
                
                length_chain.append(len(df))
            
            mean_length=np.mean(length_chain)
            
            chain_length_dict[mean_length]=chain
            
        
        longest_chain=chain_length_dict[max(chain_length_dict)]
        
    else:
        longest_chain=chains[0]
        
    
    return longest_chain
        

def get_RMSD_lovoalign(pdb1,
                       pdb2,
                       chain=None,
                       lovoalign_binary='lovoalign',
                      keep_temporary=False):
    
    """
    This function runs lovoalign and gets the RMSD value for the two aligned proteins.

    pdb1 : str
        the first pdb to be aligned
    pdb2 : str
        the second pdb to be aligned
    chain: str, default=None
        the chain of the protein that is being aligned
    lovoalign_binary : str, default=lovoalign
        used to call lovoalign--should match the command for lovoalign working computer
    keep_temporary : bool, default=False
        if True, the temporary align files will be kept rather than deleted
    """
    
    #create name for align file
    tag1=pdb1.split('/')[-1][0:4]
    tag2=pdb2.split('/')[-1][0:4]
    align_file="tmp_align_"+tag1+"_"+tag2+".log"
    
    
    # Align files using lovoalign
    if chain:
        cmd = [lovoalign_binary,"-p1",pdb1,"-c1",chain,"-p2",pdb2,"-c2",chain]

    else:
        cmd = [lovoalign_binary,"-p1",pdb1,"-p2",pdb2]
        
        
    verbose= True

    
    popen = subprocess.Popen(cmd,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT,
                         universal_newlines=True)

    f = open(align_file, 'w')
    
    for line in popen.stdout:
        f.write(line)
    f.close()

    for line in popen.stdout:
        if verbose:
            print(line,end="",flush=True)
        
    # Check for success
    return_code = popen.wait()
    if return_code != 0:
        err = "Program failed.\n"
        raise RuntimeError(err)
        
    
    #get the RMSD value from the align file
    with open(align_file) as f:
        lines = f.readlines()

    for line in lines:

        if "FINAL SCORE"  in line:
            myline=line

    line_vals=[]

    try:
        for i in myline.split(' '):

            if i != '':

                line_vals.append(i)
                
    except:
        logger.log('')

    RMSD=float(line_vals[line_vals.index('RMSD:')+1])

    #remove the align file if not keep_temporary
    if not keep_temporary:
        os.remove(align_file)
    
    
    return RMSD

def get_least_similar_pdbs (ens_dir, 
                            chains_file, 
                            out_dir,
                            copy_thermo_files=True,
                            keep_temporary=False):
    
    """
    This function finds the least similar pdbs in an ensemble and moves them into an outdirectory.

    ens_dir : str
        the directory containing the pdb files for an ensemble. These pdbs should be already synced.
    chains_file: str, default=None
        the file containing the synced chains shared between all pdbs in the ensemble.
    out_dir : str,
        the directory to where the least similar pdbs and other related files will be copied over
    copy_thermo_files : bool, default=True
        if True, ThermoMPNN ddg files for the least similar pdbs will also be copied over
    keep_temporary : bool, default=False
        if True, the temporary align files will be kept rather than deleted
    """

    #make the out_dir
    os.mkdir(out_dir)
    
    pdbs=glob.glob(ens_dir+'/*pdb')
    

    my_chain=find_longest_chain(chains_file=chains_file)
    
    RMSD_dict={}
    
    #get the RMSD values for each pair of pdbs
    for i in range(len(pdbs)):
        
        for n in range(i+1,len(pdbs)):
            
            RMSD_val=get_RMSD_lovoalign(pdbs[i],pdbs[n],chain=my_chain,keep_temporary=keep_temporary)
            RMSD_dict[float(RMSD_val)]=[pdbs[i],pdbs[n]]

    
    #find the largest RMSD value
    max_RMSD_pair = RMSD_dict[max(RMSD_dict.keys())]
    
    
    #copy the files over to the outdir
    
    for item in max_RMSD_pair:
        
        just_file_name=item.split('/')[-1]
        
        shutil.copyfile(item,out_dir+'/'+just_file_name)
        
        if copy_thermo_files == True:
            
            tag=just_file_name.split('.pdb')[0]
            
            tag_path=item.split('.pdb')[0]
            
            thermo_filepath=tag_path+'_thermo_ddg_'+my_chain+'.csv'
            
            thermo_file=tag+'_thermo_ddg_'+my_chain+'.csv'
            
            shutil.copyfile(thermo_filepath,out_dir+'/'+thermo_file)
        
    
    logger.log(max_RMSD_pair[0]+' and '+max_RMSD_pair[1]+' are the least similar pair of pdbs. They have been moved to new directory.')
    
    return max(RMSD_dict.keys())