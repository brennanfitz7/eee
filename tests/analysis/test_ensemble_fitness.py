
from eee.core import Ensemble

from eee.analysis import ensemble_fitness
from eee.core.fitness.ff import ff_on
from eee.core.fitness.ff import ff_off

import numpy as np
import pandas as pd

import pytest


def test_ensemble_fitness(ens_with_fitness):
    
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
    
    mut_energy = {"test1":0,
                  "test2":0}
    
    conditions = {"X":[0,10000],
                  "Y":[10000,0],
                  "fitness_fcn":["off","on"],
                  "select_on":"fx_obs",
                  "select_on_folded":False,
                  "temperature":1}

    # Try different fitness_fcns
    value = ensemble_fitness(ens=ens,
                             conditions=conditions,
                             mut_energy=mut_energy)
    
    assert np.array_equal(value["fitness"],[1,1])
    assert issubclass(type(value),pd.DataFrame)

    # Now select on folded. This will be folded under first condition (test2 
    # favored) and not folded under second condition (test1 favored). 
    conditions["select_on_folded"] = True
    value = ensemble_fitness(ens=ens,
                             conditions=conditions,
                             mut_energy=mut_energy)
    
    assert np.array_equal(value["fitness"],[1,0])

    conditions["fitness_fcn"] = ["off","off"]
    value = ensemble_fitness(ens=ens,
                             conditions=conditions,
                             mut_energy=mut_energy)
    
    assert np.array_equal(value["fitness"],[1,0])

    conditions["fitness_fcn"] = ["on","off"]
    value = ensemble_fitness(ens=ens,
                             conditions=conditions,
                             mut_energy=mut_energy)
    
    assert np.array_equal(value["fitness"],[0,0])

    # Check dG_obs select_on
    ens = Ensemble(gas_constant=1)
    ens.add_species(name="test1",
                    observable=True,
                    folded=False,
                    X=1)
    ens.add_species(name="test2",
                    observable=False,
                    folded=False,
                    Y=1)
    
    mut_energy = {"test1":0,
                  "test2":0}
    
    conditions = {"X":[0,1],
                  "Y":[1,0],
                  "fitness_fcn":["on","off"],
                  "select_on":"dG_obs",
                  "select_on_folded":False,
                  "temperature":1}

    value = ensemble_fitness(ens=ens,
                             conditions=conditions,
                             mut_energy=mut_energy)
    assert np.array_equal(value["fitness"],[1,2])

    
    # Now select on folded. Will be entirely unfolded --> 0,0
    conditions["select_on_folded"] = True
    value = ensemble_fitness(ens=ens,
                             conditions=conditions,
                             mut_energy=mut_energy)
    
    assert np.array_equal(value["fitness"],[0,0])


    ens = ens_with_fitness["ens"]
    conditions = ens_with_fitness["conditions"]

    value = ensemble_fitness(ens=ens,
                             conditions=conditions,
                             mut_energy={})
    assert np.array_equal(np.round(value["fitness"],2),[0.54,1.00])

    value = ensemble_fitness(ens=ens,
                             conditions=conditions,
                             mut_energy={"s1":-1.667,"s2":3.333})
    assert np.array_equal(np.round(value["fitness"],2),[0.99,0.82])
    
    value = ensemble_fitness(ens=ens,
                             conditions=conditions,
                             mut_energy={"s1":0.167,"s2":-16.667})
    assert np.array_equal(np.round(value["fitness"],2),[0.00,1.00])

    ens = Ensemble(gas_constant=1)
    ens.add_species(name="test1",
                    observable=True,
                    folded=False,
                    X=1)
    ens.add_species(name="test2",
                    observable=False,
                    folded=True,
                    Y=1)
    
    mut_energy = {"test1":0,
                  "test2":0}
    
    conditions = {"X":[0,10000],
                  "Y":[10000,0],
                  "fitness_fcn":["off","on"],
                  "select_on":"fx_obs",
                  "select_on_folded":False,
                  "temperature":1}

    value = ensemble_fitness(ens=ens,
                             conditions=conditions,
                             mut_energy=mut_energy)
    
    assert np.array_equal(value["fitness"],[1,1])
