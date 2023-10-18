"""
Read an ensemble from a file or dictionary/dataframe input.
"""

import eee
from eee.core.data import GAS_CONSTANT

import pandas as pd

import json
import os

def _search_for_key(some_dict,
                    search_key,
                    current_stack=None):
    """
    Recursively search a potentially nested dictionary for a specific key. 
    (Depth-first search.) Return the sequence of keys necessary to reach the 
    matched key. Does *not* check for duplicate keys; will return the key path
    to the first match it encounters. 
    """
    
    # Build current_stack for first iteration
    if current_stack is None:
        current_stack = []
    
    # Go through keys in dict
    for key in some_dict:
            
        # If we match the key, append to current_stack and return
        if key == search_key:
            current_stack.append(key)
            return current_stack
        
        # If we hit a nested dictionary, record the current key and search 
        # downstream for the search key. 
        if issubclass(type(some_dict[key]),dict):
            
            current_stack.append(key)
            current_stack = _search_for_key(some_dict=some_dict[key],
                                            search_key=search_key,
                                            current_stack=current_stack)
            
            # If we did not find the key, remove the key we used to access this
            # dictionary from the stack -- no match downstream. 
            if current_stack[-1] != search_key:
                current_stack = current_stack[:-1]
                
    return current_stack

def _spreadsheet_to_ensemble(df,
                             gas_constant=GAS_CONSTANT):
    """
    Load a spreadsheet and try to convert to an ensemble. Rows are treated as 
    different species; columns as keyword parameters (dG0, ligands, etc.)
    """

    #print("Loading ensemble from a spreadsheet\n",flush=True)
    
    # Get columns from spreadsheet
    df = eee.io.read_dataframe(df)
    if "name" not in df.columns:
        err = "\nspreadsheet must have a 'name' column\n\n"
        raise ValueError(err)

    # Create ensemble
    ens = eee.core.Ensemble(gas_constant=gas_constant)

    # Now load each row in as a species. 
    for idx in df.index:

        # Construct kwargs
        kwargs = {}
        for k in df.columns:
            kwargs[k] = df.loc[idx,k]
            
        # Add species
        ens.add_species(**kwargs)
    
    return ens


def _json_to_ensemble(calc_input,
                      base_path=None):
    """
    Read json and try to convert to an ensemble. 
    """

    #print("Loading ensemble from json\n",flush=True)

    if base_path is None:
        base_path = ""

    # Look for gas constant somewhere in input. If not there, use default
    key_stack = _search_for_key(calc_input,"gas_constant")
    if len(key_stack) > 0:
        gc = calc_input
        for k in key_stack:
            gc = gc[k]
        gas_constant = gc
    else:
        gas_constant = GAS_CONSTANT

    # Look for "ens" key somewhere in the json output. If it's there, pull that
    # sub-dictionary out by itself
    key_stack = _search_for_key(calc_input,"ens")
    if len(key_stack) == 0:
        err = "\njson file does not have an 'ens' key\n\n"
        raise ValueError(err)
    for k in key_stack:
        calc_input = calc_input[k]
    
    # If what remains is a string, read it as a file
    if issubclass(type(calc_input),str):
        df = os.path.join(base_path,calc_input)
        return _spreadsheet_to_ensemble(df=df,gas_constant=gas_constant)
    
    # Create ensemble from entries and validate. 
    ens = eee.core.Ensemble(gas_constant=gas_constant)
    for s in calc_input:
        
        try:
            ens.add_species(name=s,**calc_input[s])
        except TypeError:
            err = f"\nMangled json. Check ensemble keywords for species '{s}'\n\n"
            raise ValueError(err)

    return ens

def _file_to_ensemble(input_file,
                      base_path=None):
    """
    Open file and choose whether or not to read as json or spreadsheet. 
    """

    #print(f"Reading ensemble from {input_file}.",flush=True)

    # Make sure the file exists.
    if not os.path.isfile(input_file):
        err = f"\ninput_file '{input_file}' could not be read as file\n\n"
        raise FileNotFoundError(err)
    
    # If it has a .json extension
    if input_file[-5:].lower() == ".json":

        # Read json
        with open(input_file) as f: 
            input_json = json.load(f)
        
        # Construct ensemble
        ens = _json_to_ensemble(input_json,base_path=base_path)

    # Otherwise, assume it's a spreadsheet. _spreadsheet_to_ensemble uses a 
    # flexible df loader that loads files, so pass in un-edited. 
    else:
        ens = _spreadsheet_to_ensemble(df=input_file)

    return ens


def read_ensemble(input_value,
                  base_path=None):
    """
    Read an ensemble from input. The function will try to interpret the input
    as a file (json file or spreadsheet), raw json, or a pandas.DataFrame. 

    Parameters
    ----------
    input_value : 
        input to read
    base_path : str, optional
        path to any files that might be encountered when reading the input

    Notes
    -----
    If a the input is json, the ensemble is defined under the "ens" key. There
    are several acceptable formats. The simplest just lists all species, where
    the name of the species is the key and the parameters are dictionary of
    values. Any keyword arguments to the :code:`Ensemble.add_species` method 
    that are not defined revert to the default value for that argument in 
    method. Some examples follow:

    ..code-block:: json

        {
          "ens":{
            "one":{"dG0":0,"observable":true,"folded":false,"X":1},
            "two":{"dG0":0,"observable":false,"folded":false,"Y":1},
          }
        }

    You can add the special key "gas_constant" on the same level as "ens" to
    define the gas constant :

    ..code-block:: json

        {
          "gas_constant":0.00197,
          "ens":{
            "one":{"dG0":0,"observable":true,"folded":false,"X":1},
            "two":{"dG0":0,"observable":false,"folded":false,"Y":1},
          }
        }

    You can also point to a spreadsheet file defining the ensemble.

    ..code-block:: json

        {
          "gas_constant":0.00197,
          "ens":"ensemble.xlsx"
          }
        }

    You can also embed an ensemble within more complicated json defining a 
    simulation using the "ens" key. eee will search for the "ens" key and, if
    present, build the ensemble from whatever is under that key. 

    ..code-block:: json

        {
            "calc_type":"wf_tree_sim",
            "calc_params":{
                "param_1":10,
                "param_2":0.01
            },
            "gas_constant":0.00197,
            "ens":{
                "one":{"dG0":0,"observable":true,"folded":false,"X":1},
                "two":{"dG0":0,"observable":false,"folded":false,"Y":1},
            }
        }

    When reading a spreadsheet, eee treats the rows as species and the columns
    as values. It looks for columns corresponding to the keyword arguments to 
    :code:`Ensemble.add_species` and uses those values when adding each row to 
    the ensemble. Omitted keywords use their default values from the method. 
    Any columns in the spreadsheet that do not correspond to a keyword (name,
    dG0, observable, folded) are treated as ligand stoichiometries.  The
    following spreadsheet defines two species, s1 and s2. s1 interacts with
    molecule "X" with a stoichiometry of 1, s2 interacts with molecule "Y" with
    a stoichiometry of 2. 

    +------+-----+------------+---+---+
    | name | dG0 | observable | X | Y | 
    +------+-----+------------+---+---+
    | s1   | 0   | TRUE       | 1 | 0 |
    +------+-----+------------+---+---+
    | s2   | 5   | FALSE      | 0 | 2 |
    +------+-----+------------+---+---+

    A spreadsheet does NOT define a gas constant, so the default defined in 
    eee.core.data.GAS_CONSTANT is used. 
    """

    v_type = type(input_value)

    # If it's a string, parse as a file
    if issubclass(v_type,str):
        ens = _file_to_ensemble(input_value,
                                base_path=base_path)

    # If it's a dataframe, parse as a dataframe
    elif issubclass(v_type,pd.DataFrame):
        ens = _spreadsheet_to_ensemble(input_value)

    # If it's dict, treat as json
    elif issubclass(v_type,dict):
        ens = _json_to_ensemble(input_value,
                                base_path=base_path)
    
    # Otherwise, throw an error
    else: 
        err = f"\ninput_value '{input_value}' could not be read as an ensemble\n\n"
        raise ValueError(err)

    # Print status of loaded ensemble
    #print("\nBuilt the following ensemble\n")
    #print(ens.species_df)
    #print()

    return ens