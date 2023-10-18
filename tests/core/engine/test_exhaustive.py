import pytest

from eee.core.engine.exhaustive import exhaustive
from eee.core.genotype import Genotype

import numpy as np
import pandas as pd

import os
import glob

def test_exhaustive(ens_test_data,variable_types,tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    ens = ens_test_data["ens"]
    fitness_function = ens_test_data["fc"].fitness
    ddg_df = ens_test_data["ddg_df"]

    # --------------------------------------------------------------------------
    # Make sure this mutates all single sites

    gc = Genotype(ens=ens,
                  fitness_function=fitness_function,
                  ddg_df=ddg_df)


    all_muts = []    
    for site in gc.ddg_dict:
        all_muts.extend(list(gc.ddg_dict[site].keys()))
    all_muts.sort()

    df = exhaustive(gc=gc,
                    max_depth=1,
                    output_file=None,
                    return_output=True)
    
    assert len(df) == len(all_muts) + 1
    muts_seen = list(set(list(df.loc[1:,"mutations"])))
    muts_seen.sort()

    assert len(muts_seen) == len(all_muts)
    assert np.array_equal(muts_seen,all_muts)
    assert np.array_equal(df.columns,["genotype","mutations","num_mutations",
                                      "s1_ddg","s2_ddg","fitness"])


    gc = Genotype(ens=ens,
                  fitness_function=fitness_function,
                  ddg_df=ddg_df)

    df = exhaustive(gc=gc,
                    max_depth=2,
                    output_file=None,
                    return_output=True)
        
    assert np.array_equal(list(df["mutations"]),
                          ["","M1A","M1V","P2R","P2Q",
                           "M1A/P2R","M1A/P2Q","M1V/P2R","M1V/P2Q"])

    # Okay, for this one, we ask for a depth greater than we can handle (3 sites,
    # but only 2 in sites). Make sure it works fine. 

    gc = Genotype(ens=ens,
                  fitness_function=fitness_function,
                  ddg_df=ddg_df)

    df = exhaustive(gc=gc,
                    max_depth=3,
                    output_file=None,
                    return_output=True)
    
    assert np.array_equal(list(df["mutations"]),
                          ["","M1A","M1V","P2R","P2Q",
                           "M1A/P2R","M1A/P2Q","M1V/P2R","M1V/P2Q"])


    # Create a three-site system
    ddg_df = pd.DataFrame({"site":[1,1,2,2,3,3],
                           "mut":["M1A","M1V","P2R","P2Q","A3S","A3T"],
                           "s1":[1,-1,0,0,0,0],
                           "s2":[-1,1,1,0,0,0]})

    gc = Genotype(ens=ens,
                  fitness_function=fitness_function,
                  ddg_df=ddg_df)

    df = exhaustive(gc=gc,
                    max_depth=3,
                    output_file=None,
                    return_output=True)

    assert np.array_equal(list(df["mutations"]),
                          ["","M1A","M1V","P2R","P2Q","A3S","A3T",
                           "M1A/P2R","M1A/P2Q","M1V/P2R","M1V/P2Q",
                           "M1A/A3S","M1A/A3T","M1V/A3S","M1V/A3T",
                           "P2R/A3S","P2R/A3T","P2Q/A3S","P2Q/A3T",
                           "M1A/P2R/A3S","M1A/P2R/A3T",
                           "M1A/P2Q/A3S","M1A/P2Q/A3T",
                           "M1V/P2R/A3S","M1V/P2R/A3T",
                           "M1V/P2Q/A3S","M1V/P2Q/A3T"])
    

    # Create a three-site system with a single mutation at the last site
    ddg_df = pd.DataFrame({"site":[1,1,2,2,3],
                           "mut":["M1A","M1V","P2R","P2Q","A3S"],
                           "s1":[1,-1,0,0,0],
                           "s2":[-1,1,1,0,0]})

    gc = Genotype(ens=ens,
                  fitness_function=fitness_function,
                  ddg_df=ddg_df)

    df = exhaustive(gc=gc,
                    max_depth=3,
                    output_file=None,
                    return_output=True)

    assert np.array_equal(list(df["mutations"]),
                          ["","M1A","M1V","P2R","P2Q","A3S",
                           "M1A/P2R","M1A/P2Q","M1V/P2R","M1V/P2Q",
                           "M1A/A3S","M1V/A3S",
                           "P2R/A3S","P2Q/A3S",
                           "M1A/P2R/A3S",
                           "M1A/P2Q/A3S",
                           "M1V/P2R/A3S",
                           "M1V/P2Q/A3S"])

    # Now try various argument params

    ens = ens_test_data["ens"]
    fitness_function = ens_test_data["fc"].fitness
    ddg_df = ens_test_data["ddg_df"]

    gc = Genotype(ens=ens,
                  fitness_function=fitness_function,
                  ddg_df=ddg_df)
    
    # ------------------------ gc ---------------------------

    for v in variable_types["everything"]:

        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            df = exhaustive(gc=v,
                            max_depth=1,
                            output_file=None,
                            return_output=False)

    # ------------------------ max_depth ---------------------------

    # Already checked good depths. Try 0 then wacky
    df = exhaustive(gc=gc,
                    max_depth=0,
                    output_file=None,
                    return_output=True)
    assert issubclass(type(df),pd.DataFrame)
    assert len(df) == 1

    for v in variable_types["not_ints_or_coercable"]:

        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            df = exhaustive(gc=gc,
                            max_depth=v,
                            output_file=None,
                            return_output=False)
    
    with pytest.raises(ValueError):
        df = exhaustive(gc=gc,
                        max_depth=-1,
                        output_file=None,
                        return_output=False)

    # ------------------------ return_output ---------------------------
    df = exhaustive(gc=gc,
                    max_depth=1,
                    output_file=None,
                    return_output=True)
    assert issubclass(type(df),pd.DataFrame)

    df = exhaustive(gc=gc,
                    max_depth=1,
                    output_file=None,
                    return_output=False)
    assert df is None

    for v in variable_types["not_bool"]:

        try:
            if v == 0 or v == 1:
                continue
        except:
            pass

        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            df = exhaustive(gc=gc,
                            max_depth=1,
                            output_file=None,
                            return_output=v)
            
    # ------------------------ output_file ---------------------------
    df = exhaustive(gc=gc,
                    max_depth=1,
                    output_file=None,
                    return_output=True)
    assert len(glob.glob("*.csv")) == 0

    df = exhaustive(gc=gc,
                    max_depth=1,
                    output_file="junk.csv",
                    return_output=False)
    assert os.path.exists("junk.csv")


    os.chdir(current_dir)