import pytest

from eee.core import Ensemble
from eee.core.fitness.fitness import Fitness
from eee.core.fitness.ff import ff_on
from eee.core.fitness.ff import ff_off

import numpy as np
import pandas as pd

import os



def test_Fitness(): 

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

    # make sure attributes are set correctly
    fc = Fitness(ens=ens,
                 conditions=conditions)

    # fc._private_ens should be a copy of ens
    assert fc.ens is ens
    assert fc._private_ens is not ens

    assert np.array_equal(fc.temperature,[1.0,1.0])
    
    for k in fc._ligand_dict:
        assert np.array_equal(conditions[k],fc._ligand_dict[k])
    assert fc.condition_df["fitness_fcn"].iloc[0] == conditions["fitness_fcn"][0]
    assert fc.condition_df["fitness_fcn"].iloc[1] == conditions["fitness_fcn"][1]
    assert np.array_equal(fc._fitness_fcns,[ff_off,ff_on])
    assert fc._select_on == "fx_obs"
    assert np.array_equal(fc._select_on_folded,[conditions["select_on_folded"],
                                                conditions["select_on_folded"]])
    assert len(fc._fitness_kwargs) == 2
    assert issubclass(type(fc._fitness_kwargs[0]),dict)
    assert len(fc._fitness_kwargs[0]) == 0
    assert issubclass(type(fc._fitness_kwargs[1]),dict)
    assert len(fc._fitness_kwargs[0]) == 0
    assert np.array_equal(fc.temperature,np.ones(2,dtype=float))

    # We test all kinds of inputs to conditions in the test_read_fitness_conditions
    # function. This is just a pass-through that then sets attributes, so these
    # tests are sufficient. 


def test_Fitness_fitness():

    # Basic ensemble    
    ens = Ensemble(gas_constant=1)
    ens.add_species(name="test1",
                    observable=True,
                    folded=False,
                    X=1)
    ens.add_species(name="test2",
                    observable=False,
                    folded=False,
                    Y=1)
    
    conditions = {"X":[0,100],
                  "Y":[100,0],
                  "fitness_fcn":["off","on"],
                  "select_on":"fx_obs",
                  "select_on_folded":False,
                  "temperature":1}

    fc = Fitness(ens=ens,
                 conditions=conditions)

    mut_energy_array = np.array([0,0]) 
    value = fc.fitness(mut_energy_array=mut_energy_array)
    assert np.array_equal(value,[1,1])

    # Perturb so observable never populated...
    mut_energy_array = np.array([5000,0]) 
    value = fc.fitness(mut_energy_array=mut_energy_array)
    assert np.array_equal(value,[1,0])
        
    # Now set select_on_folded to True. Always unfolded.
    conditions = {"X":[0,100],
                  "Y":[100,0],
                  "fitness_fcn":["off","on"],
                  "select_on":"fx_obs",
                  "select_on_folded":True,
                  "temperature":1}

    fc = Fitness(ens=ens,
                 conditions=conditions)

    mut_energy_array = np.array([0,0]) 
    value = fc.fitness(mut_energy_array=mut_energy_array)
    assert np.array_equal(value,[0,0])

    # Perturb so observable never populated...
    mut_energy_array = np.array([5000,0]) 
    value = fc.fitness(mut_energy_array=mut_energy_array)
    assert np.array_equal(value,[0,0])


def test_Fitness_to_dict():

    # Basic ensemble    
    ens = Ensemble(gas_constant=1)
    ens.add_species(name="test1",
                    observable=True,
                    folded=False,
                    X=1)
    ens.add_species(name="test2",
                    observable=False,
                    folded=False,
                    Y=1)
    
    conditions = {"X":[0,100],
                  "Y":[100,0],
                  "fitness_fcn":["off","on"],
                  "select_on":"fx_obs",
                  "select_on_folded":False,
                  "temperature":1}

    fc = Fitness(ens=ens,conditions=conditions)

    out_dict = fc.to_dict()
    assert np.array_equal(out_dict["X"],[0,100])
    assert np.array_equal(out_dict["Y"],[100,0])
    assert np.array_equal(out_dict["fitness_fcn"],["off","on"])
    assert np.array_equal(out_dict["select_on"],["fx_obs","fx_obs"])
    assert np.array_equal(out_dict["select_on_folded"],[False,False])
    assert np.array_equal(out_dict["temperature"],[1,1])
    assert np.array_equal(out_dict["fitness_kwargs"],[{},{}])

def test_Fitness_ens():

    # Basic ensemble    
    ens = Ensemble(gas_constant=1)
    ens.add_species(name="test1",
                    observable=True,
                    folded=False,
                    X=1)
    ens.add_species(name="test2",
                    observable=False,
                    folded=False,
                    Y=1)
    
    conditions = {"X":[0,100],
                  "Y":[100,0],
                  "fitness_fcn":["off","on"],
                  "select_on":"fx_obs",
                  "select_on_folded":False,
                  "temperature":1}

    fc = Fitness(ens=ens,conditions=conditions)
    
    # fc._private_ens should be a copy of ens
    assert fc.ens is ens
    assert fc._private_ens is not ens
    
def test_Fitness_ligand_dict():

    # Basic ensemble    
    ens = Ensemble(gas_constant=1)
    ens.add_species(name="test1",
                    observable=True,
                    folded=False,
                    X=1)
    ens.add_species(name="test2",
                    observable=False,
                    folded=False,
                    Y=1)
    
    conditions = {"X":[0,100],
                  "Y":[100,0],
                  "fitness_fcn":["off","on"],
                  "select_on":"fx_obs",
                  "select_on_folded":False,
                  "temperature":1}

    fc = Fitness(ens=ens,conditions=conditions)

    assert set(fc.ligand_dict.keys()) == set(["X","Y"])
    assert np.array_equal(fc.ligand_dict["X"],[0,100])
    assert np.array_equal(fc.ligand_dict["Y"],[100,0])

def test_Fitness_select_on():

    # Basic ensemble    
    ens = Ensemble(gas_constant=1)
    ens.add_species(name="test1",
                    observable=True,
                    folded=False,
                    X=1)
    ens.add_species(name="test2",
                    observable=False,
                    folded=False,
                    Y=1)
    
    conditions = {"X":[0,100],
                  "Y":[100,0],
                  "fitness_fcn":["off","on"],
                  "select_on":"fx_obs",
                  "select_on_folded":False,
                  "temperature":1}

    fc = Fitness(ens=ens,conditions=conditions)
    assert fc.select_on == "fx_obs"

    conditions["select_on"] = "dG_obs"
    fc = Fitness(ens=ens,conditions=conditions)
    assert fc.select_on == "dG_obs"

def test_Fitness_select_on_folded():

    # Basic ensemble    
    ens = Ensemble(gas_constant=1)
    ens.add_species(name="test1",
                    observable=True,
                    folded=False,
                    X=1)
    ens.add_species(name="test2",
                    observable=False,
                    folded=False,
                    Y=1)
    
    conditions = {"X":[0,100],
                  "Y":[100,0],
                  "fitness_fcn":["off","on"],
                  "select_on":"fx_obs",
                  "select_on_folded":False,
                  "temperature":1}

    fc = Fitness(ens=ens,conditions=conditions)
    assert np.array_equal(fc.select_on_folded,[False,False])

    conditions["select_on_folded"] = True
    fc = Fitness(ens=ens,conditions=conditions)
    assert np.array_equal(fc.select_on_folded,[True,True])


def test_Fitness_fitness_kwargs():

    # Basic ensemble    
    ens = Ensemble(gas_constant=1)
    ens.add_species(name="test1",
                    observable=True,
                    folded=False,
                    X=1)
    ens.add_species(name="test2",
                    observable=False,
                    folded=False,
                    Y=1)
    
    conditions = {"X":[0,100],
                  "Y":[100,0],
                  "fitness_fcn":["off","on"],
                  "select_on":"fx_obs",
                  "select_on_folded":False,
                  "temperature":1}

    fc = Fitness(ens=ens,conditions=conditions)
    assert np.array_equal(fc.fitness_kwargs,[{},{}])

    conditions["fitness_kwargs"] = [{"rocket":1,"league":2},{}]
    fc = Fitness(ens=ens,conditions=conditions)
    assert np.array_equal(fc.fitness_kwargs,[{"rocket":1,"league":2},{}])

def test_Fitness_temperature():

    # Basic ensemble    
    ens = Ensemble(gas_constant=1)
    ens.add_species(name="test1",
                    observable=True,
                    folded=False,
                    X=1)
    ens.add_species(name="test2",
                    observable=False,
                    folded=False,
                    Y=1)
    
    conditions = {"X":[0,100],
                  "Y":[100,0],
                  "fitness_fcn":["off","on"],
                  "select_on":"fx_obs",
                  "select_on_folded":False,
                  "temperature":1}

    fc = Fitness(ens=ens,conditions=conditions)
    assert np.array_equal(fc.temperature,[1,1])

    conditions["temperature"] = np.array([10,100])
    fc = Fitness(ens=ens,conditions=conditions)
    assert np.array_equal(fc.temperature,[10,100])

def test_Fitness_condition_df():

    # Basic ensemble    
    ens = Ensemble(gas_constant=1)
    ens.add_species(name="test1",
                    observable=True,
                    folded=False,
                    X=1)
    ens.add_species(name="test2",
                    observable=False,
                    folded=False,
                    Y=1)
    
    conditions = {"X":[0,100],
                  "Y":[100,0],
                  "fitness_fcn":["off","on"],
                  "select_on":"fx_obs",
                  "select_on_folded":False,
                  "temperature":1}

    fc = Fitness(ens=ens,conditions=conditions)

    out_df = fc.condition_df
    assert issubclass(type(out_df),pd.DataFrame)
    assert len(out_df) == 2
    assert np.array_equal(out_df["fitness_fcn"],["off","on"])
    assert np.array_equal(out_df["X"],[0,100])
    assert np.array_equal(out_df["Y"],[100,0])
    assert np.array_equal(out_df["select_on"],["fx_obs","fx_obs"])
    assert np.array_equal(out_df["select_on_folded"],[False,False])
    assert np.array_equal(out_df["temperature"],[1,1])
    assert np.array_equal(out_df["fitness_kwargs"],[{},{}])

    assert len(out_df.columns) == 7