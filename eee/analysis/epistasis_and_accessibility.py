"""
Compare a deep mutational scan and epistasis calculation.
"""

import eee
from eee.calcs import read_json
from eee._private.check.compare_dict import compare_dict

import pandas as pd
import numpy as np
from tqdm.auto import tqdm

import os

def _get_rings(csv_file):
    """
    Load mutational outputs, getting all mutations seen in the first and 
    second mutational rings. 
    """
    
    print(f"Loading {csv_file}",flush=True)

    df = pd.read_csv(csv_file)
    df = df.loc[df.num_mutations<=2,:]

    ring_one_muts = []
    ring_two_genotypes = []
    ring_two_muts = []

    muts_split = df.mutations.str.split("/")
    for i in tqdm(range(1,len(muts_split))):
        muts_split[i].sort()
        ring_two_muts.extend(muts_split[i])
        muts_split[i] = "/".join(muts_split[i])

    df.loc[:,"mutations"] = muts_split

    ring_one_muts = set(df.loc[df.num_mutations == 1,"mutations"])
    ring_two_genotypes = set(df.loc[df.num_mutations == 2,"mutations"])
    ring_two_muts = set(ring_two_muts)
        
    return df, ring_one_muts, ring_two_muts, ring_two_genotypes
        
        
def epistasis_and_accessibility(dms_dir,
                                accessible_dir):
    """
    Compare a deep-mutational scan to an accessibility calculation. Identify 
    double mutants that have a fitness greater than or less than the wildtype
    that are inaccessible due to passing through a non-fit genotype. Return the
    epistasis in fitness for these genotypes. 
    
    Parameters
    ----------
    dms_dir : str
        directory holding a dms calculation (must have simulation.json)
    accessible_dir : str
        directory holding an accessible_paths calculation (must have 
        simulation.json) 
    
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

    # Read acc_calc and make sure it is correct calc type
    acc_calc, acc_kwargs = read_json(os.path.join(accessible_dir,"input","simulation.json"))
    if acc_calc.calc_type != "accessible_paths":
        err = f"\naccessible_dir calc_type is '{acc_calc.calc_type}', not 'accessible_paths'\n\n"
        raise ValueError(err)

    # Make sure the two calculations have the same ensembles
    if not compare_dict(dms_calc.ens.to_dict(),
                        acc_calc.ens.to_dict()):
        err = "\ndms and accessible ensembles are not identical\n\n"
        raise ValueError(err)

    # Make sure the two calculations have the same fitness conditions
    if not compare_dict(dms_calc.fc.to_dict(),
                        acc_calc.fc.to_dict()):
        err = "\ndms and accessible fitness conditions are not identical\n\n"
        raise ValueError(err)
    
    # Make sure both calcs have at least two mutation max_depth
    if dms_kwargs["max_depth"] < 2:
        err = "\ndms max_depth must be 2 or higher to do this analysis\n\n"
        raise ValueError(err)

    if acc_kwargs["max_depth"] < 2:
        err = "\naccesibility max_depth must be 2 or higher to do this analysis\n\n"
        raise ValueError(err)
    
    # Get all genotypes with one and two mutations (from dms)
    dms_file = os.path.join(dms_dir,dms_kwargs["output_file"])
    dms_df, dms_ring_one_muts, dms_ring_two_muts, dms_ring_two_genotypes = _get_rings(dms_file)

    # Get accessible genotypes with one and two mutation (from acc)
    acc_file = os.path.join(accessible_dir,acc_kwargs["output_file"])
    acc_df, acc_ring_one_muts, acc_ring_two_muts, acc_ring_two_genotypes = _get_rings(acc_file)
                    
    # Get fitness values for all single mutants as a dictionary
    dms_ring_one = dms_df.loc[dms_df.num_mutations == 1,:]
    single_fitness = dict(zip(dms_ring_one.mutations,
                            dms_ring_one.fitness))

    num_dms_single = len(dms_ring_one_muts)
    num_dms_double = len(dms_ring_two_genotypes)

    num_acc_single = len(acc_ring_one_muts)
    num_acc_double = len(acc_ring_two_genotypes)

    num_acc_in_two_not_one = len(acc_ring_two_muts - acc_ring_one_muts)

    # Get wildtype fitness
    f00 = np.array(dms_df.loc[dms_df.num_mutations == 0,"fitness"])[0]

    # Get double mutants that are inaccessible, but have a fitness >= wildtype
    dms_ok_two_mask = np.logical_and(dms_df.fitness >= f00,
                                     dms_df.num_mutations == 2)
    dms_ok_two = set(dms_df.loc[dms_ok_two_mask,"mutations"])
    
    inacc = set(dms_ok_two) - set(acc_ring_two_genotypes)
    
    inacc_mask = dms_df.mutations.isin(inacc)

    # Create dictionary to hold outputs
    ep_out = ["m01","m10","m11",
              "f00","f01","f10","f11",
              "mag","sign01","sign10","ep_class",
              "stuck"]
    ep_out = dict([(e,[]) for e in ep_out])

    # Go through inaccessible mutants
    print("Measuring epistasis for inaccesible mutants.",flush=True)
    for idx in tqdm(dms_df.index[inacc_mask]):

        # Get data
        row = dms_df.loc[idx,:]
            
        # Get single mutants that make this up
        singles = row.mutations.split("/")
        
        # Get their fitness values
        f01 = single_fitness[singles[0]]
        f10 = single_fitness[singles[1]]
        f11 = row.fitness
        
        # Calculate the epistasis in fitness between these mutations
        mag, sign01, sign10, ep_class = eee.analysis.epistasis.get_epistasis(f00,
                                                                             f01,
                                                                             f10,
                                                                             f11)
        
        # Record the epistasis
        ep_out["m01"].append(singles[0])
        ep_out["m10"].append(singles[1])
        ep_out["m11"].append(row.mutations)
        ep_out["f00"].append(f00)
        ep_out["f01"].append(f01)
        ep_out["f10"].append(f10)
        ep_out["f11"].append(f11)
        ep_out["mag"].append(mag)
        ep_out["sign01"].append(sign01)
        ep_out["sign10"].append(sign10)
        ep_out["ep_class"].append(ep_class)
        
        if ((f11 < f10) and (f11 < f01)) or \
           ((f11 < f10) and (f01 < f00)) or \
           ((f11 < f01) and (f10 < f00)):
            ep_out["stuck"].append(True)
        else:
            ep_out["stuck"].append(False)

    # Final dataframe
    ep_df = pd.DataFrame(ep_out)

    single_pct = num_acc_single/num_dms_single*100
    print(f"{num_acc_single} of {num_dms_single} single mutants ({single_pct:.4f}%) are accessible")

    double_pct = num_acc_double/num_dms_double*100
    print(f"{num_acc_double} of {num_dms_double} double mutants ({double_pct:.4f}%) are accessible")

    print(f"{num_acc_in_two_not_one} mutations are accessible in ring two but not ring one")
    print(f"{len(inacc)} double mutants are not accessible despite having fitness >= wildtype")
    print()

    num_stuck = np.sum(ep_df["stuck"])
    num_valley = np.sum(np.logical_not(ep_df["stuck"]))
    print(f"+ {num_stuck} get stuck on a higher-fitness single mutant")
    print(f"+ {num_valley} would have to cross a fitness valley")
    print("",flush=True)

    return ep_df