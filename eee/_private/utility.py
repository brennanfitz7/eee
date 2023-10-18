"""
Functions for miscellaneous tasks. 
"""

import numpy as np

import copy

def prep_for_json(some_dict,work_on_copy=True):
    """
    Prepare a dictionary for writing to json. Converts any np.ndarray to lists
    and any numpy datatypes to standard types. 

    Paramters
    ---------
    some_dict : dict
        dictionary to write

    Returns
    -------
    cleaned_dict : dict
        cleaned up dictionary
    """

    if not issubclass(type(some_dict),dict):
        err = f"\nsome_dict {some_dict} should be a dictionary\n\n"
        raise ValueError(err)

    if work_on_copy:
        some_dict = copy.deepcopy(some_dict)

    # Go through a dictionary...
    for k, v in some_dict.items():     

        # If the value is a dictionary...
        v_type = type(v)   
        if issubclass(v_type, dict):
            prep_for_json(v,work_on_copy=False)

        # Otherwise...
        else:
            
            # If this is a non-dict, non-string iterable
            if not issubclass(v_type,str) and hasattr(v,"__iter__"):

                # Coerce numpy array to list
                if issubclass(v_type,np.ndarray):
                    v = list(v)

                # Go through each element in the list and coerce to a
                # built in type
                for i in range(len(v)):
                    local_v_type = type(v[i])
                    if np.issubdtype(local_v_type,np.integer):
                        v[i] = int(v[i])
                    elif np.issubdtype(local_v_type,np.floating):
                        v[i] = float(v[i])
                    elif np.issubdtype(local_v_type,np.bool_):
                        v[i] = bool(v[i])
                    else:
                        continue

            # If not an iterable, coerce individual types
            else:
                if np.issubdtype(v_type,np.integer):
                    v = int(v)
                elif np.issubdtype(v_type,np.floating):
                    v = float(v)
                elif np.issubdtype(v_type,np.bool_):
                    v = bool(v)
                else:
                    continue

            some_dict[k] = v

    return some_dict