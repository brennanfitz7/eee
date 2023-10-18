"""
Load and validate eee json files.
"""

from eee.io.read_ddg import read_ddg
from eee.io.read_tree import read_tree
from eee.io.read_ensemble import read_ensemble
from eee.calcs import CALC_AVAILABLE

import json
import inspect
import os

def _validate_calc_kwargs(calc_type,
                          calc_function,
                          kwargs):
    """
    Make sure the kwargs defined in the json file match the arguments for the
    calc_function. This does not check the types of the arguments (that is done
    within each class) but does make sure we have the correct argument names.
    It will generate a human-readable error message if the arguments are 
    not correct.
    """

    # Set of kwargs found
    kwargs_found = set(list(kwargs.keys()))
        
    # Get signature with all kwargs
    sig = inspect.signature(calc_function)

    # Grab arguments, separating them into those that must be defined 
    # and those with defaults. 
    required = []
    have_defaults = []
    for param in sig.parameters:
        if sig.parameters[param].default is inspect._empty:
            if param != "self":
                required.append(param)
        else:
            have_defaults.append(param)

    # Create sets of arguments in signature
    required = set(required)
    have_defaults = set(have_defaults)
    all_allowed_args = required | have_defaults

    # Look for args that are not found or extra args
    missing_required = required - required.intersection(kwargs)
    extra_args = kwargs_found - kwargs_found.intersection(all_allowed_args)
    
    # Optimistically start by assuming success
    success = True
    miss_err = ""
    extra_err = ""
    
    # Missing arguments. Create a human-readable error string
    if len(missing_required) > 0:
        success = False

        missing_required = list(missing_required)
        missing_required.sort()

        miss_err = "\nThe following required keys are not defined:\n"
        for m in missing_required:
            miss_err += f"    {m}\n"
        miss_err += "\n"

    # Extra arguments.. Create a human-readable error string
    if len(extra_args) > 0:
        success = False

        extra_args = list(extra_args)
        extra_args.sort()

        extra_err = "\nThe following keys are not allowed:\n"
        for e in extra_args:
            extra_err += f"    {e}\n"
        extra_err += "\n"

    # If we failed above, construct a human-readable error string
    if not success:

        # Start with main error
        err = "\nThe json file does not have the correct arguments for calc_type\n"
        err += f"'{calc_type}'.\n\n"
        err = err + miss_err + extra_err 
        
        name = f"{calc_function}"
        dashes = len(name)*"-"
        
        err += f"\ncalc_type '{calc_type}' details:\n\n{name}\n{dashes}\n"   
        
        # Drop whole doc string into error message
        err += f"{calc_function.__doc__}\n\n"

        raise ValueError(err)
    
    return kwargs

def read_json(json_file,
              use_stored_seed=False):
    """
    Load a json file describing a simulation. This file must have the 
    following top-level keys:

        "calc_type" : a string indicating what kind of calculation this is. 
        "ens" : a dictionary of species describing the ensemble.

    The other allowed top-level keys are:

        "conditions" : a dictionary describing conditions at which to do 
                       fitness calculations
        "ddg_df" : spreadsheet file with the effects of mutations on each 
                   conformation in the ensemble. 
        "calc_params" : any parameters needed to run the calculation indicated 
                        by calc_type.
        "seed" : an integer describing the seed for reproducible calcs. 
        "gas_constant" : a positive float holding the gas constant (sets the 
                         energy units for the calculation)

    Parameters
    ----------
    json_file : str
        json file to load
    use_stored_seed : bool, default=False
        The 'seed' key in the json file (if present) is ignored unless
        use_stored_seed is set to True. The only time to re-use the seed 
        is to restart a simulation or reproduce it exactly for testing 
        purposes. 

    Returns
    -------
    sc : SimulationContainer subclass
        initialized SimulationContainer subclass with ensemble, fitness, and
        ddg loaded.
    calc_params : dict
        dictionary with run parameters. sc.run(**calc_params) will run the 
        calculation. 
    """

    # Read json file
    with open(json_file) as f:
        calc_input = json.load(f)

    base_path = os.path.dirname(os.path.abspath(json_file))

    setup_kwargs = {}
    
    # ---------------------------- calc_type -----------------------------------
    if "calc_type" not in calc_input:
        err = "\njson must have 'calc_type' key in top level that defines the\n"
        err += "calculation being done.\n\n"
        raise ValueError(err)
    calc_type = calc_input.pop("calc_type")

    if not issubclass(type(calc_type),str) or calc_type not in CALC_AVAILABLE:
    
        err = f"\ncalc_type '{calc_type}' is not recognized. calc_type should\n"
        err += "be one of:\n"
        for a in CALC_AVAILABLE:
            err += f"    {a}\n"
        err += "\n"
        raise ValueError(err)

    calc_class = CALC_AVAILABLE[calc_type]
    
    # ---------------------------- calc_params ---------------------------------

    # Make sure we have a "calc_params" key, even if empty
    if "calc_params" not in calc_input:
        calc_input["calc_params"] = {}

    # Load newick here so we don't have to keep track of the file when/if we
    # start a simulation in new directory
    if "tree" in calc_input["calc_params"]:
        newick_file = os.path.join(base_path,calc_input["calc_params"]["tree"])
        calc_input["calc_params"]["tree"] = read_tree(newick_file)

    # -------------------------------- ens -------------------------------------

    # Create an ensemble from the 'ens' key.
    if "ens" not in calc_input:
        err = "\njson must have 'ens' that defines the thermodynamic ensemble.\n\n"
        raise ValueError(err)
    
    # Read the ensemble
    ens = read_ensemble(input_value=calc_input,
                        base_path=base_path)
    
    setup_kwargs["ens"] = ens

    # ---------------------------- conditions ----------------------------------

    if "conditions" in calc_input:
        
        conditions = calc_input["conditions"]
        if issubclass(type(conditions),str):
            conditions_file = os.path.join(base_path,conditions)
            if not os.path.isfile(conditions_file):
                err = "\nconditions entry '{conditions}' could not be read\n\n"
                raise ValueError(err)
            conditions = conditions_file

        setup_kwargs["conditions"] = conditions
    
    # ------------------------------ ddg_df ------------------------------------
    # Load ddg_df here so we don't have to keep track of the file when/if we
    # start a simulation in new directory
    if "ddg_df" in calc_input:
        ddg_file = os.path.join(base_path,calc_input["ddg_df"])
        ddg_df = read_ddg(ddg_file)
        setup_kwargs["ddg_df"] = ddg_df
    
    # ------------------------------- seed -------------------------------------

    # Drop the seed unless we are requesting it to be kept. 
    if use_stored_seed and "seed" in calc_input:
        setup_kwargs["seed"] = calc_input["seed"]

    # ------------------------------- seed -------------------------------------

    calc_params = _validate_calc_kwargs(calc_type=calc_type,
                                        calc_function=calc_class.__init__,
                                        kwargs=setup_kwargs)

    # Set up the calculation class. 
    sc = calc_class(**setup_kwargs)

    calc_params = _validate_calc_kwargs(calc_type=calc_type,
                                        calc_function=sc.run,
                                        kwargs=calc_input["calc_params"])

    return sc, calc_params


    

