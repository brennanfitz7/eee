"""
Load an rcsb file (pdb or cif) into a pandas data frame, pre-processing to 
remove solvent atoms, etc. 
"""

from eee.core.data import AA_3TO1

import pandas as pd
import numpy as np

def _read_structure_pdb(pdb_file):
    """
    Load the coordinates from a pdb file into a pandas data frame. Should 
    generally be called via read_structure.
    """
    
    out = {"model":[],
           "class":[],
           "chain":[],
           "resid":[],
           "resid_num":[],
           "alternate":[],
           "elem":[],
           "atom":[],
           "atom_num":[],
           "x":[],
           "y":[],
           "z":[],
           "b":[],
           "occ":[]}

    model = None
    with open(pdb_file) as f:
        for line in f:
            if line[0:6] in ["ATOM  ","HETATM"]:
                
                if model is None:
                    model = 1

                try:
                    out["model"].append(model)
                    out["class"].append(line[0:6].strip())
                    out["chain"].append(line[21])
                    out["atom_num"].append(int(line[6:11]))
                    out["atom"].append(line[12:16].strip())
                    out["resid"].append(line[17:20].strip())
                    out["resid_num"].append(line[22:28].strip())

                    alternate = line[16]
                    if alternate == " ":
                        alternate = "."

                    out["alternate"].append(alternate)
                    out["x"].append(float(line[30:38]))
                    out["y"].append(float(line[38:46]))
                    out["z"].append(float(line[46:54]))
                    out["occ"].append(float(line[54:60]))
                    out["b"].append(float(line[60:66]))
                    out["elem"].append(line[75:].strip())

                except Exception as e:
                    print(f"Could not parse line:\n\n{line}\n\n",flush=True)
                    raise(e)
                
            if line.startswith("MODEL"):
                model = int(line[6:].strip())



    return pd.DataFrame(out)


def _read_structure_cif(cif_file):
    """
    Load the coordinates from a cif file into a pandas dataframe. Should 
    generally be called via read_structure.
    """

    out = {"model":[],
           "class":[],
           "chain":[],
           "resid":[],
           "resid_num":[],
           "alternate":[],
           "elem":[],
           "atom":[],
           "atom_num":[],
           "x":[],
           "y":[],
           "z":[],
           "b":[],
           "occ":[]}

    with open(cif_file) as f:
        for line in f:

            columns = line.split()
            if columns[0] in ["ATOM","HETATM"]:
                
                try:
                    out["model"].append(int(columns[20]))
                    out["class"].append(columns[0].strip())
                    out["atom_num"].append(int(columns[1]))
                    out["atom"].append(columns[3])
                    out["alternate"].append(columns[4].strip())
                    out["resid"].append(columns[5].strip())
                    out["chain"].append(columns[6].strip())
                    out["resid_num"].append(columns[8])
                    out["x"].append(float(columns[10]))
                    out["y"].append(float(columns[11]))
                    out["z"].append(float(columns[12]))
                    out["occ"].append(float(columns[13]))
                    out["b"].append(float(columns[14]))
                    out["elem"].append(columns[2].strip())
                except Exception as e:
                    print(f"Could not parse line:\n\n{line}\n\n",flush=True)
                    raise(e)

    return pd.DataFrame(out)
    

def read_structure(rcsb_file,
                   remove_solvent=True,
                   remove_non_protein_polymer=True,
                   remove_multiple_models=True,
                   remove_alternate_conf=True,
                   remove_hydrogens=True):
    """
    Load an rcsb file (pdb or cif) into a pandas data frame.
    
    Parameters
    ----------
    rcsb_file : str
        rcsb file with .pdb or .cif extension
    remove_solvent : bool, default=True
        remove HOH and SOL atoms
    remove_non_protein_polymer : bool, default=True
        remove non-protein polymers
    remove_multiple_models : bool, default=True
        remove multiple models if present, taking only first
    remove_alternate_conf : bool, default=True
        remove alternate conformations, taking only "A" if present
    remove_hydrogens : bool, default=True
        remove hydrogen atoms
    
    Returns
    -------
    df : pandas.DataFrame
        pandas dataframe with coordinates and other atom information contained
        in the rcsb file. 
    """
    
    # Figure out which file type this is and load
    if rcsb_file[-4:] == ".cif":
        df = _read_structure_cif(rcsb_file)
    elif rcsb_file[-4:] == ".pdb":
        df = _read_structure_pdb(rcsb_file)
    else:
        err = "file type not recognized. Should be .cif or .pdb\n"
        raise ValueError(err)

    # Remove waters
    if remove_solvent:
        df = df.loc[np.logical_not(np.logical_and(df["class"] == "HETATM",
                                                  df["resid"].isin(["HOH","SOL"]))),:]
    # Remove non-protein polymers
    if remove_non_protein_polymer:
        aa_resids = list(AA_3TO1.keys())
        
        non_protein_mask = np.logical_and(df["class"] == "ATOM",
                                          np.logical_not(df["resid"].isin(aa_resids)))
        keep_mask = np.logical_not(non_protein_mask)
        df = df.loc[keep_mask,:]
        
    # Remove all models except first
    if remove_multiple_models:
        if len(np.unique(df.model)) > 1:
            df = df.loc[df.model == df.model.iloc[0],:]
            
    # Take first alternate conformation at each residue
    if remove_alternate_conf:
        df = df.loc[df.alternate.isin([".","A"]),:]
        
    # Remove hydrogens
    if remove_hydrogens:
        df = df.loc[df["elem"] != "H",:]

    # Clean up dataframe indexes
    df.index = np.arange(len(df.index),dtype=int)

    return df
        