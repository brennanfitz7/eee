
from eee.io.read_structure import _read_structure_cif
from eee.io.read_structure import _read_structure_pdb
from eee.io.read_structure import read_structure

import pandas as pd
import numpy as np

def test__read_structure_cif(test_cifs):
    
    for k in test_cifs:
        df = _read_structure_cif(test_cifs[k])
        assert issubclass(type(df),pd.DataFrame)

def test__read_structure_pdb(test_pdbs):
    
    for k in test_pdbs:
        df = _read_structure_pdb(test_pdbs[k])
        assert issubclass(type(df),pd.DataFrame)

def test_read_structure(test_cifs,test_pdbs):

    # Head-to-head comparison of pdb and cif loading
    for k in test_cifs:

        # Load from cif and pdb. All except atom_num, resid_num and chain should
        # be identical
        cif_df = read_structure(test_cifs[k])
        cif_df = cif_df.drop(columns=["atom_num","resid_num","chain"])
        pdb_df = read_structure(test_pdbs[f"{k[:-4]}.pdb"])
        pdb_df = pdb_df.drop(columns=["atom_num","resid_num","chain"])

        # Should be same length
        assert len(cif_df) == len(pdb_df)

        # Should be identical
        print("comparing cif and pdb for",k)
        assert np.array_equal(np.array(cif_df),np.array(pdb_df))

    # Make sure remove_solvent flag is working like we think
    df = read_structure(test_cifs["2bb2.cif"],remove_solvent=True)
    assert len(df.loc[df["resid"] == "HOH",:]) == 0

    df = read_structure(test_cifs["2bb2.cif"],remove_solvent=False)
    assert len(df.loc[df["resid"] == "HOH",:]) > 0

    # Make sure remove_non_protein_polymer flag is working like we think
    df = read_structure(test_cifs["1lcc.cif"],remove_non_protein_polymer=True)
    assert len(df.loc[df["resid"] == "DA",:]) == 0

    df = read_structure(test_cifs["1lcc.cif"],remove_non_protein_polymer=False)
    assert len(df.loc[df["resid"] == "DA",:]) > 0

    # make sure remove_multiple_models flag is working like we think
    df = read_structure(test_cifs["7ui5.cif"],remove_multiple_models=True)
    assert len(np.unique(df["model"])) == 1

    df = read_structure(test_cifs["7ui5.cif"],remove_multiple_models=False)
    assert len(np.unique(df["model"])) == 10

    # Make sure remove_alternate_conf flag is working like we think. 
    df = read_structure(test_cifs["5qu4.cif"],remve_alternate_conf=True)
    alt = np.unique(df["alternate"])
    alt.sort()
    assert np.array_equal(alt,[".","A"])

    # Make sure remove_alternate_conf flag is working like we think. 
    df = read_structure(test_cifs["5qu4.cif"],remove_alternate_conf=False)
    alt = np.unique(df["alternate"])
    alt.sort()
    assert np.array_equal(alt,[".","A","B"])

    # Make sure remove_hydrogens flag is working like we think
    df = read_structure(test_cifs["7ui5.cif"],remove_hydrogens=True)
    assert len(df.loc[df["elem"] == "H",:]) == 0

    df = read_structure(test_cifs["7ui5.cif"],remove_hydrogens=False)
    assert len(df.loc[df["elem"] == "H",:]) > 0