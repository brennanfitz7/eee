
from eee.io.read_structure import read_structure
from eee.io.write_pdb import write_pdb

from eee._private.interface import create_new_dir
from eee._private.interface import launch
from eee._private.interface import rmtree

import pandas as pd
import numpy as np

import os
import shutil

def clean_structure(df,
                    foldx_binary="foldx",
                    verbose=False,
                    keep_temporary=False,
                    remove_multiple_models=True):
    """
    Run a structure through foldx to build missing atoms in sidechains. This 
    will delete residues with incomplete backbones.

    Parameters
    ----------
    df : pandas.DataFrame
        dataframe loaded from an rscb structure file
    foldx_binary : str
        foldx_binary to use
    verbose : bool, default=False
        write out all output to standard output
    keep_temporary : bool, default=False
        do not delete temporary files
    remove_multiple_models : bool
        bool that determines whether multiple models are kept in the pdb file after read_structure

    Returns
    -------
    df : pandas.DataFrame
        dataframe with residues cleaned up by foldx.
    """

    tmp_dir = create_new_dir()    
    write_pdb(df,os.path.join(tmp_dir,"input.pdb"))

    cmd = [foldx_binary,
            "-c","PDBFile",
            "--fixSideChains","1",
            "--pdb","input.pdb"]
    
    launch(cmd=cmd,
           run_directory=tmp_dir,
           suppress_output=(not verbose))
    
    shutil.move(os.path.join(tmp_dir,"PF_input.fxout"),
                os.path.join(tmp_dir,"tmp-output.pdb"))
    
    new_df = read_structure(os.path.join(tmp_dir,"tmp-output.pdb"),remove_multiple_models=remove_multiple_models)

    print('this is new df after cleaning',new_df)

    # foldx will drop all hetatms. bring them back in
    hetatm_df = df.loc[df["class"] == "HETATM",:]
    new_df = pd.concat((new_df,hetatm_df))
    new_df.index = np.arange(len(new_df.index),dtype=int)
    
    if not keep_temporary:
        rmtree(tmp_dir)

    return new_df