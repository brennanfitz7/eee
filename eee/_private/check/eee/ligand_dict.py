"""
Check the sanity of ligand_dict. Returns a validated ligand_dict and the length of
the ligand_dict conditions. 
"""

from eee._private.check.standard import check_float

import numpy as np
import pandas as pd

def check_ligand_dict(ligand_dict):
    """
    Check the sanity of ligand_dict, making sure all values are arrays with the 
    same length.
    
    Parameters
    ----------
    ligand_dict : dict
        dictionary holding chemical potentials keyed to species name

    Returns
    -------
    ligand_dict : dict
        validated ligand_dict, where all values are arrays of the same length
    length : int
        the length of the ligand_dict conditions
    """

    # Should be a dictionary
    if not issubclass(type(ligand_dict),dict) or issubclass(type(ligand_dict),type):
        err = "ligand_dict should be a dictionary that keys chemical species to chemical potential\n"
        raise ValueError(err)
    
    # Empty dict: allowed, just return
    if len(ligand_dict) == 0:
        return ligand_dict, 1

    # Check each value...
    ligand_lengths = []
    for lig in ligand_dict:

        # Make sure not disallowed value class that has __iter__
        value_type = type(ligand_dict[lig])
        for bad in [type,pd.DataFrame,dict]:
            if issubclass(value_type,bad):
                err = f"\nligand_dict['{lig}'] cannot be type {bad}\n\n"
                raise ValueError(err)
        
        # If it is iterable, we have work to do to check types
        if hasattr(ligand_dict[lig],"__iter__"):

            # Make sure there is actually something in the iterable
            if len(ligand_dict[lig]) == 0:
                err = f"\nligand_dict['{lig}'] must have a length > 0\n\n"
                raise ValueError(err)

            # If a string, try to coerce into a float.
            if issubclass(value_type,str):
                v = check_float(ligand_dict[lig],variable_name=f"ligand_dict['{lig}']")
                ligand_dict[lig] = np.ones(1,dtype=float)*v
            
            # If we get here, coerce into a numpy array.
            else:
                ligand_dict[lig] = np.array(ligand_dict[lig],dtype=float)
                if len(ligand_dict[lig]) != 1:
                    ligand_lengths.append(len(ligand_dict[lig]))

        else:

            # Single value float
            v = check_float(value=ligand_dict[lig],
                            variable_name=f"ligand_dict['{lig}']") 
            ligand_dict[lig] = np.ones(1,dtype=float)*v

    # All lengths are 1: return
    if len(ligand_lengths) == 0:
        return ligand_dict, 1

    # Unique non-one ligand_lengths
    ligand_lengths = list(set(ligand_lengths))
    if len(ligand_lengths) > 1:
        err = "\nall values in ligand_dict must have the same length\n\n"
        raise ValueError(err)
    
    # Take any values with length one and make them the same length as the 
    # longer value. 
    final_ligand_length = ligand_lengths[0]
    for lig in ligand_dict:
        if len(ligand_dict[lig]) == 1:
            ligand_dict[lig] = np.ones(final_ligand_length,dtype=float)*ligand_dict[lig][0]

    return ligand_dict, final_ligand_length