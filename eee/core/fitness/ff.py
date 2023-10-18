"""
Fitness functions for a specific observable. Any function in this file that 
starts with ff_ will be made available for evolutionary simulations. These 
functions must take the value of an observable as their first argument. They 
may also take further keyword arguments (passed in via fitness_kwargs). Whatever
occurs after ff_ will be the name of the function in things like json files 
using the function. 

.. code-block::python

    def ff_example_function(value,*args,**kwargs):
        '''
        Return the fitness of a value as that value times multiply_by. 
        '''

        multiply_by = kwargs["multiply_by]

        return value*multiply_by

We could use this in a simulation set up as follows. This set of kwargs will
apply the function ff_example_function to "fx_obs." ff_example_function will 
multiply "fx_obs" by 5 at an iptg chemical potential of 1. 

.. code-block::python

    {"fitness_fcns":["example_function"],
    "select_on": "fx_obs",
    "fitness_kwargs":{"multiply_by":5},
    "ligand_dict":{"iptg":1}}
    
"""

def ff_on(value,*args,**kwargs):
    """
    Fitness is linearly proportional to value. Useful for simulating selection
    to keep observable 'on'.
    """
    return value

def ff_off(value,*args,**kwargs):
    """
    Fitness is linearly proportional to 1 - value. Useful for simulating
    selection to keep observable 'off'. 
    """
    return 1 - value

def ff_neutral(value,*args,**kwargs):
    """
    Fitness is always 1.0, modeling no selection on observable. 
    """
    return 1.0

def ff_on_above(value,*args,**kwargs):
    """
    Fitness is 1.0 when observable above a threshold, 0.0 when below it. When
    used, `fitness_kwargs = {"threshold":VALUE}` must be set. 
    """

    threshold = kwargs["threshold"]
    if value >= threshold:
        return 1.0
    
    return 0.0

def ff_on_below(value,*args,**kwargs):
    """
    Fitness is 0.0 when observable above a threshold, 1.0 when below it. When
    used, `fitness_kwargs = {"threshold":VALUE}` must be set. 
    """

    threshold = kwargs["threshold"]
    if value <= threshold:
        return 1.0
    
    return 0.0


def _get_ff_available():
    
    ff_available = {}
    
    possible = dict(globals())
    for p in possible:
        this_poss = possible[p]
        if callable(this_poss):
            if p.startswith("ff_"):
                ff_available[p[3:]] = this_poss

    return ff_available
                
# Register calculations available
FF_AVAILABLE = _get_ff_available()