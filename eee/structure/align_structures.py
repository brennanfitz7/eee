
from eee.io.write_pdb import write_pdb
from eee.io.read_structure import read_structure
from eee._private import logger

from eee._private.interface import create_new_dir
from eee._private.interface import launch
from eee._private.interface import rmtree

import numpy as np

import os

def align_structures(dfs,
                     lovoalign_binary="lovoalign",
                     verbose=False,
                     keep_temporary=False):
    """
    Align the structures in the dataframes using lovoalign. Align all structures
    individually to the first dataframe.
    
    Parameters
    ----------
    dfs : list
        list of pandas dataframes containing structures
    lovoalign_binary : default = "lovoalign"
        lovoalign binary
    verbose : bool, default = True
        whether or not to print lovoalign output
    keep_temporary : bool, default=False
        do not delete temporary files
    
    Returns
    -------
    dfs : list
        list of pandas dataframes containing structures with x, y, z coordinates
        updated
    """
    
    # No structures to align
    if len(dfs) < 2:
        return dfs
    
    tmp_dir = create_new_dir()

    # Write out pdb files to align
    files = []
    
    for i, df in enumerate(dfs):
        no_het_df = df.loc[df["class"] == "ATOM",:]

        files.append(f"tmp-align_{i}.pdb")
        write_pdb(no_het_df,os.path.join(tmp_dir,files[-1]))
         
    # Go through all but first file
    for i, f in enumerate(files[1:]):
        
        try:
            # Align file to first file using lovoalign
            out_file = f"tmp-out_{i+1}.pdb"
            cmd = [lovoalign_binary,"-p1",f,"-p2",files[0],"-o",out_file]
            launch(cmd=cmd,run_directory=tmp_dir)
                    
            # Move coordinates from aligned structure into dfs
            new_df = read_structure(os.path.join(tmp_dir,out_file))
            mask = dfs[i+1]["class"] == "ATOM"

            dfs[i+1].loc[mask,"x"] = np.array(new_df["x"])
            dfs[i+1].loc[mask,"y"] = np.array(new_df["y"])
            dfs[i+1].loc[mask,"z"] = np.array(new_df["z"])
        except:
            logger.log('A file failed at lovoalign alignment.')
        
    # Remove all temporary files
    if not keep_temporary:
        rmtree(tmp_dir)
        
    return dfs
