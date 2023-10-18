"""
Functions for mapping fitness functions in ff.py to string names. 
"""

from .ff import FF_AVAILABLE

def map_fitness_fcn(value,return_as):
    """
    Take an input value and map back and forth between the name of that function
    and a pointer to the function.
    
    Parameters
    ----------
    value : function or str
        should be a fitness function or string name of a function
    return_as : str
        must be either "function" (return the function) or "string" (return 
        the string name of the function). 

    Returns
    -------
    result : function or str
        function or str from mapping, depending on return_as
    """
    
    str_to_fcn = FF_AVAILABLE
    fcn_to_str = dict([(str_to_fcn[s],s) for s in str_to_fcn])

    # If what came in is a function ...
    if callable(value):
        fcn = value
        if value in fcn_to_str:
            fcn_str = fcn_to_str[fcn]
        else:
            fcn_str = f"{fcn}"

    # If what came in is a string ...
    elif issubclass(type(value),str):
        fcn_str = value.strip().lower()
        if value in str_to_fcn:
            fcn = str_to_fcn[fcn_str]
        else:
            fcn = None

    # If neither, die
    else:
        err = f"\nvalue '{value}' should be a fitness function or string\n\n"
        raise ValueError(err)

    # Check return_as argument
    if not issubclass(type(return_as),str) or issubclass(type(return_as),type):
        err = f"\nreturn_as ({return_as}) should be 'function' or 'string'\n\n"
        raise ValueError(err)

    # Return function if possible
    if return_as == "function":

        if fcn is None:
            err = f"fitness function string ('{fcn_str}') should be one of:\n"
            for k in str_to_fcn:
                err += f"    {k}\n"
            err += "\n"
            raise ValueError(err)
        
        return fcn
    
    # Return string. (If a function was passed in that is not in
    # FF_AVAILABLE, just return it's name. Allows a custom function
    # to be passed in). 
    elif return_as == "string":
        return fcn_str
    
    # Or die. 
    else:
        err = f"\nreturn_as ({return_as}) should be 'function' or 'string'\n\n"
        raise ValueError(err)
