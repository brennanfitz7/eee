import pytest

from eee.core import Ensemble
from eee.core.fitness.read_fitness_conditions import read_fitness_conditions

import numpy as np
import pandas as pd

import os

def test_read_fitness_conditions(tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    # -------------------------------------------------------------------------
    # Basic ensemble 

    ens = Ensemble(gas_constant=1)
    ens.add_species(name="test1",
                    observable=True,
                    folded=False,
                    X=1)
    ens.add_species(name="test2",
                    observable=False,
                    folded=True,
                    Y=1)
    
    conditions = {"X":[0,10000],
                  "Y":[10000,0],
                  "fitness_fcn":["off","on"],
                  "select_on":"fx_obs",
                  "select_on_folded":True,
                  "temperature":1}
    
    df, ligand_dict = read_fitness_conditions(ens=ens,
                                      conditions=conditions,
                                      default_fitness_kwargs=None,
                                      default_select_on="fx_obs",
                                      default_select_on_folded=True,
                                      default_temperature=298.15)

    assert np.array_equal(df.loc[:,"X"],conditions["X"])
    assert np.array_equal(df.loc[:,"Y"],conditions["Y"])
    assert np.array_equal(df.loc[:,"fitness_fcn"],conditions["fitness_fcn"])
    assert np.array_equal(df.loc[:,"select_on"],["fx_obs","fx_obs"])
    assert np.array_equal(df.loc[:,"select_on_folded"],[True,True])
    assert np.array_equal(df.loc[:,"temperature"],[1,1])
    assert np.array_equal(df.loc[:,"fitness_kwargs"],[{},{}])
    
    assert np.array_equal(ligand_dict["X"],conditions["X"])
    assert np.array_equal(ligand_dict["Y"],conditions["Y"])
    assert len(ligand_dict) == 2

    # Pass in {} for default fitness_kwargs
    df, ligand_dict = read_fitness_conditions(ens=ens,
                                      conditions=conditions,
                                      default_fitness_kwargs={},
                                      default_select_on="fx_obs",
                                      default_select_on_folded=True,
                                      default_temperature=298.15)

    # -------------------------------------------------------------------------
    # Ensemble has reserved name for a ligand

    ens = Ensemble(gas_constant=1)
    ens.add_species(name="test1",
                    observable=True,
                    folded=False,
                    select_on=1)
    
    conditions = {"X":[0,10000],
                  "Y":[10000,0],
                  "fitness_fcn":["off","on"],
                  "select_on":"fx_obs",
                  "select_on_folded":True,
                  "temperature":1}
    
    # Send in ligand with a reserved name
    with pytest.raises(ValueError):
        df, ligand_dict = read_fitness_conditions(ens=ens,
                                            conditions=conditions,
                                            default_fitness_kwargs=None,
                                            default_select_on="fx_obs",
                                            default_select_on_folded=True,
                                            default_temperature=298.15)

    # -------------------------------------------------------------------------
    # Missing required column

    conditions = {"X":[0,10000],
                  "Y":[10000,0],
                  "select_on":"fx_obs",
                  "select_on_folded":True,
                  "temperature":1}
    
    with pytest.raises(ValueError):
        df, ligand_dict = read_fitness_conditions(ens=ens,
                                          conditions=conditions,
                                          default_fitness_kwargs=None,
                                          default_select_on="fx_obs",
                                          default_select_on_folded=True,
                                          default_temperature=298.15)




    # -------------------------------------------------------------------------
    # pass in conditions with non-ligand columns

    ens = Ensemble(gas_constant=1)
    ens.add_species(name="test1",
                    observable=True,
                    folded=False,
                    X=1)
    ens.add_species(name="test2",
                    observable=False,
                    folded=True,
                    Y=1)
    
    conditions = {"X":[0,10000],
                  "Y":[10000,0],
                  "not_really_there":[1,2],
                  "fitness_fcn":["off","on"],
                  "select_on":"fx_obs",
                  "select_on_folded":True,
                  "temperature":1}
    with pytest.raises(ValueError):
        df, ligand_dict = read_fitness_conditions(ens=ens,
                                          conditions=conditions,
                                          default_fitness_kwargs=None,
                                          default_select_on="fx_obs",
                                          default_select_on_folded=True,
                                          default_temperature=298.15)

    # -------------------------------------------------------------------------
    #  make sure columns get default values. 

    ens = Ensemble(gas_constant=1)
    ens.add_species(name="test1",
                    observable=True,
                    folded=False,
                    X=1)
    ens.add_species(name="test2",
                    observable=False,
                    folded=True,
                    Y=1)
    
    conditions = {"X":[0,10000],
                  "Y":[10000,0],
                  "fitness_fcn":["off","on"]}
    
    df, ligand_dict = read_fitness_conditions(ens=ens,
                                      conditions=conditions,
                                      default_fitness_kwargs=None,
                                      default_select_on="dG_obs",
                                      default_select_on_folded=False,
                                      default_temperature=300)

    assert np.array_equal(df.loc[:,"X"],conditions["X"])
    assert np.array_equal(df.loc[:,"Y"],conditions["Y"])
    assert np.array_equal(df.loc[:,"fitness_fcn"],conditions["fitness_fcn"])
    assert np.array_equal(df.loc[:,"select_on"],["dG_obs","dG_obs"])
    assert np.array_equal(df.loc[:,"select_on_folded"],[False,False])
    assert np.array_equal(df.loc[:,"temperature"],[300,300])
    assert np.array_equal(df.loc[:,"fitness_kwargs"],[{},{}])
    
    assert np.array_equal(ligand_dict["X"],conditions["X"])
    assert np.array_equal(ligand_dict["Y"],conditions["Y"])
    assert len(ligand_dict) == 2

    # -------------------------------------------------------------------------
    # Make sure ligands get assigned if not in conditions

    ens = Ensemble(gas_constant=1)
    ens.add_species(name="test1",
                    observable=True,
                    folded=False,
                    X=1)
    ens.add_species(name="test2",
                    observable=False,
                    folded=True,
                    Y=1)
    
    conditions = {"fitness_fcn":["off","on"],
                  "select_on":"fx_obs",
                  "select_on_folded":True,
                  "temperature":1}
    
    df, ligand_dict = read_fitness_conditions(ens=ens,
                                      conditions=conditions,
                                      default_fitness_kwargs=None,
                                      default_select_on="fx_obs",
                                      default_select_on_folded=True,
                                      default_temperature=298.15)

    assert np.array_equal(df.loc[:,"X"],[0,0])
    assert np.array_equal(df.loc[:,"Y"],[0,0])
    assert np.array_equal(df.loc[:,"fitness_fcn"],conditions["fitness_fcn"])
    assert np.array_equal(df.loc[:,"select_on"],["fx_obs","fx_obs"])
    assert np.array_equal(df.loc[:,"select_on_folded"],[True,True])
    assert np.array_equal(df.loc[:,"temperature"],[1,1])
    assert np.array_equal(df.loc[:,"fitness_kwargs"],[{},{}])
    
    assert np.array_equal(ligand_dict["X"],[0,0])
    assert np.array_equal(ligand_dict["Y"],[0,0])
    assert len(ligand_dict) == 2

    # -------------------------------------------------------------------------
    # Test parsing of fitness_fcn

    ens = Ensemble(gas_constant=1)
    ens.add_species(name="test1",
                    observable=True,
                    folded=False,
                    X=1)
    ens.add_species(name="test2",
                    observable=False,
                    folded=True,
                    Y=1)
    
    conditions = {"X":[0,1,2],
                  "Y":[0,1,2],
                  "fitness_fcn":["off","on","neutral"]}
    
    df, ligand_dict = read_fitness_conditions(ens=ens,
                                      conditions=conditions)

    assert np.array_equal(df.loc[:,"X"],conditions["X"])
    assert np.array_equal(df.loc[:,"Y"],conditions["Y"])
    assert np.array_equal(df.loc[:,"fitness_fcn"],conditions["fitness_fcn"])
    assert np.array_equal(df.loc[:,"fitness_kwargs"],[{},{},{}])

    # not real function
    conditions["fitness_fcn"] = ["off","on","not_real"]
    with pytest.raises(ValueError):
        read_fitness_conditions(ens=ens,conditions=conditions)

    # real function, not enough kwargs
    conditions["fitness_fcn"] = ["off","on","on_above"]
    with pytest.raises(ValueError):
        read_fitness_conditions(ens=ens,conditions=conditions)
    
    # Test the test, should work now
    conditions["fitness_fcn"] = ["off","on","off"]
    read_fitness_conditions(ens=ens,conditions=conditions)

    # -------------------------------------------------------------------------
    # Test parsing of fitness_kwargs

    ens = Ensemble(gas_constant=1)
    ens.add_species(name="test1",
                    observable=True,
                    folded=False,
                    X=1)
    ens.add_species(name="test2",
                    observable=False,
                    folded=True,
                    Y=1)
    
    conditions = {"X":[0,1,2],
                  "fitness_fcn":["off","on","neutral"],
                  "fitness_kwargs":[{},{},{}]}
    
    df, ligand_dict = read_fitness_conditions(ens=ens,conditions=conditions)
    assert np.array_equal(df.loc[:,"fitness_kwargs"],[{},{},{}])

    # now pass in default fitness_kwargs as {}
    conditions["fitness_kwargs"] = {}
    df, ligand_dict = read_fitness_conditions(ens=ens,conditions=conditions)
    assert np.array_equal(df.loc[:,"fitness_kwargs"],[{},{},{}])

    conditions["fitness_kwargs"] = {"threshold":1}
    df, ligand_dict = read_fitness_conditions(ens=ens,conditions=conditions)
    assert np.array_equal(df.loc[:,"fitness_kwargs"],[{"threshold":1},
                                                      {"threshold":1},
                                                      {"threshold":1}])

    # Send in useless fitness_kwargs
    conditions = {"X":[0,1],
                  "fitness_fcn":["on","on_above"],
                  "fitness_kwargs":[{},{}]}
    with pytest.raises(ValueError):
        df, ligand_dict = read_fitness_conditions(ens=ens,conditions=conditions)        

    # Send in unparsable fitness_kwargs
    conditions = {"X":[0,1],
                  "fitness_fcn":["on","on"],
                  "fitness_kwargs":[4,5]}
    with pytest.raises(ValueError):
        df, ligand_dict = read_fitness_conditions(ens=ens,conditions=conditions)        


    # Send in usable fitness_kwargs
    conditions = {"X":[0,1],
                  "fitness_fcn":["on","on_above"],
                  "fitness_kwargs":[{},{"threshold":0.1}]}
    df, ligand_dict = read_fitness_conditions(ens=ens,conditions=conditions)
    assert np.array_equal(df.loc[:,"fitness_kwargs"],[{},
                                                      {"threshold":0.1}])  

    # Send in usable fitness_kwargs ("on" ignores threshold)
    conditions = {"X":[0,1],
                  "fitness_fcn":["on","on_above"],
                  "fitness_kwargs":{"threshold":0.1}}
    df, ligand_dict = read_fitness_conditions(ens=ens,conditions=conditions)
    assert np.array_equal(df.loc[:,"fitness_kwargs"],[{"threshold":0.1},
                                                      {"threshold":0.1}])  

    # -------------------------------------------------------------------------
    # Test parsing of select_on

    ens = Ensemble(gas_constant=1)
    ens.add_species(name="test1",
                    observable=True,
                    folded=False,
                    X=1)
    ens.add_species(name="test2",
                    observable=False,
                    folded=True,
                    Y=1)
    
    conditions = {"X":[0,10000],
                  "fitness_fcn":["off","on"],
                  "select_on":"fx_obs"}
    df, ligand_dict = read_fitness_conditions(ens=ens,conditions=conditions)
    assert np.array_equal(df["select_on"],["fx_obs","fx_obs"])

    conditions["select_on"] = "dG_obs"
    df, ligand_dict = read_fitness_conditions(ens=ens,conditions=conditions)
    assert np.array_equal(df["select_on"],["dG_obs","dG_obs"])

    conditions["select_on"] = ["dG_obs","dG_obs"]
    df, ligand_dict = read_fitness_conditions(ens=ens,conditions=conditions)
    assert np.array_equal(df["select_on"],["dG_obs","dG_obs"])

    # should fail because must be the same across all conditions
    conditions["select_on"] = ["dG_obs","fx_obs"]
    with pytest.raises(ValueError):
        df, ligand_dict = read_fitness_conditions(ens=ens,conditions=conditions)
    
    # Should fail because select_on must be a string
    conditions["select_on"] = 1.0
    with pytest.raises(ValueError):
        df, ligand_dict = read_fitness_conditions(ens=ens,conditions=conditions)

    # Should fail because select_on must be an interpretable string
    conditions["select_on"] = "not_recognized_string"
    with pytest.raises(ValueError):
        df, ligand_dict = read_fitness_conditions(ens=ens,conditions=conditions)

    # -------------------------------------------------------------------------
    # Test parsing of select_on_folded

    ens = Ensemble(gas_constant=1)
    ens.add_species(name="test1",
                    observable=True,
                    folded=False,
                    X=1)
    ens.add_species(name="test2",
                    observable=False,
                    folded=True,
                    Y=1)
    
    conditions = {"X":[0,10000],
                  "fitness_fcn":["off","on"],
                  "select_on_folded":True}
    df, ligand_dict = read_fitness_conditions(ens=ens,conditions=conditions)
    assert np.array_equal(df["select_on_folded"],[True,True])

    conditions["select_on_folded"] = [False,True]
    df, ligand_dict = read_fitness_conditions(ens=ens,conditions=conditions)
    assert np.array_equal(df["select_on_folded"],[False,True])

    # Should work -- interpretable as bool
    conditions["select_on_folded"] = [0,1]
    df, ligand_dict = read_fitness_conditions(ens=ens,conditions=conditions)
    assert np.array_equal(df["select_on_folded"],[False,True])

    # should not work, not interpetable as bool
    conditions["select_on_folded"] = "not_a_bool"
    with pytest.raises(ValueError):
        df, ligand_dict = read_fitness_conditions(ens=ens,conditions=conditions)

    # -------------------------------------------------------------------------
    # Test parsing of temperature

    ens = Ensemble(gas_constant=1)
    ens.add_species(name="test1",
                    observable=True,
                    folded=False,
                    X=1)
    ens.add_species(name="test2",
                    observable=False,
                    folded=True,
                    Y=1)
    
    conditions = {"X":[0,10000],
                  "fitness_fcn":["off","on"],
                  "temperature":100}
    df, ligand_dict = read_fitness_conditions(ens=ens,conditions=conditions)
    assert np.array_equal(df["temperature"],[100,100])

    conditions["temperature"] = [1,100]
    df, ligand_dict = read_fitness_conditions(ens=ens,conditions=conditions)
    assert np.array_equal(df["temperature"],[1,100])

    # should not work, not interpretable as float
    conditions["temperature"] = "not_a_float"
    with pytest.raises(ValueError):
        df, ligand_dict = read_fitness_conditions(ens=ens,conditions=conditions)

    # should not work, bad float
    conditions["temperature"] = -1.5
    with pytest.raises(ValueError):
        df, ligand_dict = read_fitness_conditions(ens=ens,conditions=conditions)

    # -------------------------------------------------------------------------
    # Test ligand_dict construction

    ens = Ensemble(gas_constant=1)
    ens.add_species(name="test1",
                    observable=True,
                    folded=False,
                    X=1)
    ens.add_species(name="test2",
                    observable=False,
                    folded=True,
                    Y=1)
    
    conditions = {"X":[0,10000],
                  "fitness_fcn":["off","on"]}
    df, ligand_dict = read_fitness_conditions(ens=ens,conditions=conditions)
    assert np.array_equal(ligand_dict["X"],[0,10000])
    assert np.array_equal(ligand_dict["Y"],[0,0])

    # Should work. Can be coerced to float
    conditions = {"X":[0,"10000"],
                  "fitness_fcn":["off","on"]}
    df, ligand_dict = read_fitness_conditions(ens=ens,conditions=conditions)
    assert np.array_equal(ligand_dict["X"],[0,10000])
    assert np.array_equal(ligand_dict["Y"],[0,0])

    # Should not work. Cannot be coerced to float
    conditions = {"X":[0,"test"],
                  "fitness_fcn":["off","on"]}
    with pytest.raises(ValueError):
        df, ligand_dict = read_fitness_conditions(ens=ens,conditions=conditions)

        
    # -------------------------------------------------------------------------
    # read conditions from file
   
    ens = Ensemble(gas_constant=1)
    ens.add_species(name="test1",
                    observable=True,
                    folded=False,
                    X=1)
    ens.add_species(name="test2",
                    observable=False,
                    folded=True,
                    Y=1)
    
    conditions = {"X":[0,10000],
                  "Y":[10000,0],
                  "fitness_fcn":["off","on"],
                  "select_on":"fx_obs",
                  "select_on_folded":True,
                  "temperature":1}
    input_df = pd.DataFrame(conditions)
    input_df.to_csv("test.csv")

    df, ligand_dict = read_fitness_conditions(ens=ens,conditions="test.csv")
    
    assert np.array_equal(df.loc[:,"X"],conditions["X"])
    assert np.array_equal(df.loc[:,"Y"],conditions["Y"])
    assert np.array_equal(df.loc[:,"fitness_fcn"],conditions["fitness_fcn"])
    assert np.array_equal(df.loc[:,"select_on"],["fx_obs","fx_obs"])
    assert np.array_equal(df.loc[:,"select_on_folded"],[True,True])
    assert np.array_equal(df.loc[:,"temperature"],[1,1])
    assert np.array_equal(df.loc[:,"fitness_kwargs"],[{},{}])
    
    assert np.array_equal(ligand_dict["X"],conditions["X"])
    assert np.array_equal(ligand_dict["Y"],conditions["Y"])
    assert len(ligand_dict) == 2

    # -------------------------------------------------------------------------
    # read conditions from dataframe
   
    ens = Ensemble(gas_constant=1)
    ens.add_species(name="test1",
                    observable=True,
                    folded=False,
                    X=1)
    ens.add_species(name="test2",
                    observable=False,
                    folded=True,
                    Y=1)
    
    conditions = {"X":[0,10000],
                  "Y":[10000,0],
                  "fitness_fcn":["off","on"],
                  "select_on":"fx_obs",
                  "select_on_folded":True,
                  "temperature":1}
    input_df = pd.DataFrame(conditions)

    df, ligand_dict = read_fitness_conditions(ens=ens,conditions=input_df)
    
    assert np.array_equal(df.loc[:,"X"],conditions["X"])
    assert np.array_equal(df.loc[:,"Y"],conditions["Y"])
    assert np.array_equal(df.loc[:,"fitness_fcn"],conditions["fitness_fcn"])
    assert np.array_equal(df.loc[:,"select_on"],["fx_obs","fx_obs"])
    assert np.array_equal(df.loc[:,"select_on_folded"],[True,True])
    assert np.array_equal(df.loc[:,"temperature"],[1,1])
    assert np.array_equal(df.loc[:,"fitness_kwargs"],[{},{}])
    
    assert np.array_equal(ligand_dict["X"],conditions["X"])
    assert np.array_equal(ligand_dict["Y"],conditions["Y"])
    assert len(ligand_dict) == 2

    os.chdir(current_dir)