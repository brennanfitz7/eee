from eee.io import read_structure
from eee.data import AA_3TO1

import numpy as np
import os

def write_fasta(df,fasta_file:str,overwrite=False):

    """
    Write a pdb file given a pandas dataframe.

    Parameters
    ----------
    df : pandas.DataFrame
        dataframe with structural data (generally created using read_structure)
    fasta_file : str
        name of fasta file to write
    overwrite : bool, default=False
        overwrite the pdb file if it exists
    """

    if os.path.exists(fasta_file):
        if not overwrite:
            err = f"fasta_file {fasta_file} already exists.\n"
            raise FileExistsError(err)
        else:
            if os.path.isfile(fasta_file):
                os.remove(fasta_file)
            else:
                err = f"fasta_file {fasta_file} exists but is not a regular file.\n"
                err += "Cannot overwrite.\n"
                raise FileExistsError(err)

    seq_name=fasta_file.split('.')[0].split('/')[-1] ##remove splits later? may be redundant?

    mask = np.logical_and(df.atom == "CA",
                          df["class"] == "ATOM")
    this_df = df.loc[mask,:]        
    seq = [AA_3TO1[aa] for aa in this_df["resid"]]

    # Write fasta
    with open(seq_name+'.fasta','w') as f:
        f.write(f">{seq_name}\n{''.join(seq)}\n")
