"""
Compare a deep mutational scan and epistasis calculation.
"""

import eee
from eee.calcs import read_json
from eee._private.check.standard import check_bool

import pandas as pd
import numpy as np
from tqdm.auto import tqdm

import os

def dms_epistasis(dms_dir,
                  allow_neutral=True):
    """
    Analyze epistasis and path accessibility from a deep mutational scan.
    
    Parameters
    ----------
    dms_dir : str
        directory holding a dms calculation (must have simulation.json)
    allow_neutral : bool, default=True
        treat mutations that have no effect on fitness as accessible
    
    Returns
    -------
    ep_df : pandas.DataFrame
        dataframe holding epistasis for mutations that are not accessible but 
        have fitness >= wildtype
    """
    
    # Read dms_calc and make sure it is correct calc type
    dms_calc, dms_kwargs = read_json(os.path.join(dms_dir,"input","simulation.json"))
    if dms_calc.calc_type != "dms":
        err = f"\ndms_dir calc_type is '{dms_calc.calc_type}', not 'dms'\n\n"
        raise ValueError(err)
    
    # Make sure calculation had at least two mutations
    if dms_kwargs["max_depth"] < 2:
        err = "\ndms max_depth must be 2 or higher to do this analysis\n\n"
        raise ValueError(err)
    
    dms_file = os.path.join(dms_dir,dms_kwargs["output_file"])

    # Decide how to treat allowable steps
    allow_neutral = check_bool(allow_neutral,
                                variable_name="allow_neutral")
    if allow_neutral:
        operator = np.greater_equal
    else:
        operator = np.greater

    # Load dms into dataframe
    dms_df = pd.read_csv(dms_file)
    dms_df = dms_df.loc[dms_df.num_mutations<=2,:]
                
    # Get fitness values for all single mutants as a dictionary
    dms_ring_one = dms_df.loc[dms_df.num_mutations == 1,:]
    single_fitness = dict(zip(dms_ring_one.mutations,
                              dms_ring_one.fitness))

    # Get wildtype fitness
    f00 = np.array(dms_df.loc[dms_df.num_mutations == 0,"fitness"])[0]

    # Get only doubles
    dms_df = dms_df.loc[dms_df["num_mutations"]==2,:]
    num_doubles = len(dms_df.index)
    
    # Build output dictionary
    ep_out = {}

    str_keys = ["m10","m01","m11","ep_class","cycle_type"]
    for k in str_keys:
        ep_out[k] = ["" for _ in range(num_doubles)]

    float_keys = ["f00","f10","f01","f11","mag"]
    float_template = np.zeros(num_doubles,dtype=float)
    for k in float_keys:
        ep_out[k] = float_template.copy()

    bool_keys = ["sign10","sign01","accessible","allowable"]
    bool_template = np.zeros(num_doubles,dtype=bool)
    for k in bool_keys:
        ep_out[k] = bool_template.copy()

    int_keys = ["num_valleys","num_barriers"]
    int_template = np.zeros(num_doubles,dtype=int)
    for k in int_keys:
        ep_out[k] = int_template.copy()

    pbar = tqdm(total=num_doubles)

    with pbar:

        for i, idx in enumerate(dms_df.index):

            m11 = dms_df.loc[idx,"mutations"]
            m10, m01 = m11.split("/")
            
            f10 = single_fitness[m10]
            f01 = single_fitness[m01]
            f11 = dms_df.loc[idx,"fitness"]

            mag, sign10, sign01, ep_class = eee.analysis.epistasis.get_epistasis(f00,
                                                                                 f10,
                                                                                 f01,
                                                                                 f11)
            # Legs
            leg_00_10 = operator(f10,f00)
            leg_10_11 = operator(f11,f10)
            leg_00_01 = operator(f01,f00)
            leg_01_11 = operator(f11,f01)

            # Cycle as string of 1s and 0s
            cycle_type = [int(leg_00_01),
                          int(leg_00_10),
                          int(leg_10_11),
                          int(leg_01_11)]
            cycle_type = "c" + "".join([f"{v}" for v in cycle_type])

            # Can we get across space by one of the paths?
            accessible = np.logical_or(np.logical_and(leg_00_10,leg_10_11),
                                       np.logical_and(leg_00_01,leg_01_11))

            # If we could get to the double, would it be fit enough relative to 
            # wildtype?
            allowable = operator(f11,f00)

            # Valleys are drops in fitness going from wt to single
            num_valleys = (1 - int(leg_00_01)) + (1 - int(leg_00_10))

            # Barriers are drops in fitness going single to double
            num_barriers = (1 - int(leg_10_11)) + (1 - int(leg_01_11))

            # Record the epistasis
            ep_out["m10"][i] = m10
            ep_out["m01"][i] = m01
            ep_out["m11"][i] = m11
            ep_out["f00"][i] = f00
            ep_out["f01"][i] = f01
            ep_out["f10"][i] = f10
            ep_out["f11"][i] = f11
            ep_out["mag"][i] = mag
            ep_out["sign10"][i] = sign10
            ep_out["sign01"][i] = sign01
            ep_out["ep_class"][i] = ep_class
            ep_out["cycle_type"][i] = cycle_type
            ep_out["accessible"][i] = accessible
            ep_out["allowable"][i] = allowable
            ep_out["num_valleys"][i] = num_valleys
            ep_out["num_barriers"][i] = num_barriers
        
            pbar.update()

    # Create final dataframe with columns in useful order
    ep_out = pd.DataFrame(ep_out)

    column_order = ["m10","m01","m11",
                    "f00","f10","f01","f11",
                    "mag","sign10","sign01","ep_class",
                    "allowable","accessible",
                    "num_valleys","num_barriers","cycle_type"]

    ep_out = ep_out.loc[:,column_order]

    return ep_out
