"""
Validate temperature inputs.
"""

from eee._private.check.standard import check_float

import numpy as np

def check_temperature(temperature,num_conditions):
    """
    Check the sanity of a temperature, returning a numpy array num_conditions 
    long holding temperature. 

    Parameters
    ----------
    temperature : float or numpy.ndarray
        float with temperature; if numpy array, must be num_conditions long
    num_conditions : int
        number of conditions for calculation. will either expand T to this 
        length or must match the length of temperature

    Returns
    -------
    temperature : numpy.ndarray
        array of floats num_conditions long
    """

    if issubclass(type(temperature),type):
        err = f"\n temperature ({temperature}) cannot be a type\n\n"
        raise ValueError(err)

    if hasattr(temperature,"__iter__") and not issubclass(type(temperature),str):

        if issubclass(type(temperature),dict):
            err = "\nT cannot be a dictionary\n\n"
            raise ValueError(err)

        if len(temperature) != num_conditions:
            err = "\ntemperature should be the same length as the number of conditions. "
            err += f"{len(temperature)} vs. {num_conditions})\n\n"
            raise ValueError(err)

        # Convert to floats
        temperature = list(temperature)
        for i in range(len(temperature)):
            temperature[i] = check_float(temperature[i],
                                         variable_name=f"temperature[{i}]",
                                         minimum_allowed=0,
                                         minimum_inclusive=False)

        # Coerce to numpy array
        temperature = np.array(temperature,dtype=float)
    
    else:
        temperature = check_float(temperature,
                                  variable_name="temperature",
                                  minimum_allowed=0,
                                  minimum_inclusive=False)
        temperature = temperature*np.ones(num_conditions,dtype=float)

    return temperature
