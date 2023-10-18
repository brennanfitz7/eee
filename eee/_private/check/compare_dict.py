
import numpy as np
import pandas as pd

def compare_dict(dict_one,dict_two):
    """
    Compare two dictionaries, returning True if they are identical and False
    if they differ in any key or value. Handles dict, set, list-like, numpy,
    single-value, and pandas DataFrame comparisons.

    Parameters
    ----------
    dict_one : dict
        first dictionary
    dict_two: dict
        second dictionary
        
    Returns
    -------
    identical : bool
        True or False if two dictionaries are equivalent. 
    """
    
    # dict_one not a dict
    if not issubclass(type(dict_one),dict):
        return False

    # dict_two not a dict
    if not issubclass(type(dict_two),dict):
        return False
    
    # Keys are not identical between dicts
    if set(dict_one.keys()) != set(dict_two.keys()):
        return False
        
    # Go through keys in dict_one
    for key in dict_one:
                
        # Comparing dictionaries
        if issubclass(type(dict_one[key]),dict):
            
            # dict_one[key] is a dict but dict_two[key] is not
            if not issubclass(type(dict_two[key]),dict):
                return False
                  
            # Recursively compare two dicts
            status = compare_dict(dict_one[key],
                                  dict_two[key])
            
            # If comparison fails, return False
            if not status:
                return False
        
        # Comparing sets
        elif issubclass(type(dict_one[key]),set):
            
            if not issubclass(type(dict_two[key]),set):
                return False
            
            if dict_one[key] != dict_two[key]:
                return False
            
        # Comparing pandas dataframe
        elif issubclass(type(dict_one[key]),pd.DataFrame):

            if not issubclass(type(dict_two[key]),pd.DataFrame):
                return False
        
            # Make sure columns are identical
            if len(dict_one[key].columns) != len(dict_two[key].columns):
                return False
            if not np.array_equal(dict_one[key].columns,
                                  dict_two[key].columns):
                return False
        
            # Make sure content is identical
            if not np.array_equal(dict_one[key],dict_two[key]):
                return False
            
        # Comparing iterable but not dict, set of DataFrame
        elif hasattr(dict_one[key],"__iter__"):
            
            # Not iterable
            if not hasattr(dict_two[key],"__iter__"):
                return False
            
            # Check equivalence
            if not np.array_equal(dict_one[key],dict_two[key]):
                return False
            
        # Comparing any other data type
        else:

            # Nulls (which also catches NA) are not the same
            if pd.isnull(dict_one[key]):
                return False
            
            if pd.isnull(dict_two[key]):
                return False
        
            if dict_one[key] != dict_two[key]:
                return False
            
    # If we get here, everything is the same. Return True. 
    return True
            