from eee.io import read_dataframe

import pandas as pd
import numpy as np

import copy

def check_dataframe(value,variable_name=None):
    """
    Validate a dataframe, converting dataframe-like things into pandas
    dataframes, and checking for errors along the way. 

    Parameters
    ----------
    value : 
        object to convert into a dataframe. See Notes for details.
    variable_name : str, optional   
        name of variable to help with user error messages

    Returns
    -------
    df : pandas.DataFrame
        pandas dataframe representation of the input

    Notes
    -----
    This function can take a variety of input types. If the input is a 
    pandas.DataFrame, this is returned (by reference; not a copy). If the input
    is a string, the function attempts to read this in as a file. If the input
    is a dictionary, the function passes it as direct input to the pd.DataFrame
    init function. If the input is list-like, each element must be a dictionary.
    The elements of the list are treated as rows; the dictionary keys are 
    columns; the dictionary values are the value to place in that column/row
    coordinate. Each column must either be seen in all rows OR in exactly one 
    row. If a column is only seen in one row, its value is filled for all other 
    rows. 

    Example:

    The input :code:`[{"A":1,"B":2},{"B":3}]` would become:
    
    +---+---+
    | A | B |
    +===+===+
    | 1 | 2 | 
    +---+---+
    | 1 | 3 |
    +---+---+ 
    """

    # Create helpful text for error message
    if variable_name is not None:
        err_base = f"{variable_name} {value}"
    else:
        err_base = f"{value}"
    
    # Check for types, which cause extreme wackiness (edge case)
    if issubclass(type(value),type):
        err = f"\n{err_base} cannot be a type\n\n"
        raise ValueError(err)
    
    # If already a dataframe, return a dataframe
    if issubclass(type(value),pd.DataFrame):
        return value

    # If a string, try to load as a file
    if issubclass(type(value),str):
        return read_dataframe(df=value,
                              remove_extra_index=True)

    # If a dictionary, try to convert this into a dataframe directly
    if issubclass(type(value),dict):

        value_is_iter = [k for k in value if hasattr(value[k],"__iter__")]
        value_is_dict = [k for k in value if issubclass(type(value[k]),dict)]

        # If true, this is a mixture of dict and other iterable. Wrap the dicts
        # in lists so dict are treated as a singleton values instead of rows of
        # values 
        if len(value_is_iter) != len(value_is_dict):

            # Get the number of rows by looking at length of one of the non-dict
            # entries
            non_dict = list(set(value_is_iter) - set(value_is_dict))
            num_rows = len(value[non_dict[0]])

            # work on a copy because we are modifying before passing into 
            # pd.DataFrame
            value = copy.deepcopy(value)
            for k in value:

                # Replace any dict with list of dicts
                if issubclass(type(value[k]),dict):
                    value[k] = [value[k] for _ in range(num_rows)]

        return pd.DataFrame(value)

    # Must be an iterable at this point
    if not hasattr(value,"__iter__"):
        err = f"\n{err_base} must be a list of dicts or a dict\n"
        raise ValueError(err)


    # Build a column dictionary. 
    all_columns = {}
    for i, row in enumerate(value):

        if not issubclass(type(row),dict):
            err = f"\n{err_base} must be a list of dicts or a dict\n"
            raise ValueError(err)
        
        for key in row.keys():
            if key not in all_columns:
                all_columns[key] = [None for _ in range(len(value))]

            all_columns[key][i] = row[key]

    # Expand columns as needed.
    for c in all_columns:
        
        num_none = np.sum([v is None for v in all_columns[c]])
        if num_none == 0:
            continue
        elif num_none == len(value) - 1:
            v = [v for v in all_columns[c] if v is not None][0]
            all_columns[c] = [v for _ in range(len(value))]
        else:
            err = f"\n{err_base} has mis-matched numbers of values\n\n"
            raise ValueError(err)
        
    return pd.DataFrame(all_columns)
        
