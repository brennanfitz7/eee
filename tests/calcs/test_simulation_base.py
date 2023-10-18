import pytest

from eee.calcs.simulation_base import Simulation
from eee.core import Fitness
from eee.core import Genotype
from eee.io.read_ensemble import read_ensemble
from eee.core.fitness.read_fitness_conditions import read_fitness_conditions

import numpy as np
import pandas as pd

import os
import json

class Simulation_no_run(Simulation):
    """
    Bad dummy class. Does not have run defined, so it should throw 
    NotImplementedError.
    """

    calc_type = "fake_no_run"

class Simulation_no_calc_type(Simulation):
    """
    Bad dummy class. Does not have calc_type defined, so it should throw 
    NotImplementedError.
    """

    def run(self,*args,**kwargs):
        pass

class SimulationTester(Simulation):
    """
    Simulation cannot be run on its own. This dummy subclass can 
    be used to test its core functionality.
    """

    calc_type = "fake"

    def run(self,
            output_directory,
            population_size=100,
            mutation_rate=0.01,
            num_generations=100,
            write_prefix="eee_sim",
            write_frequency=1000):
        
        calc_params = {"population_size":population_size,
                       "mutation_rate":mutation_rate,
                       "num_generations":num_generations,
                       "write_prefix":write_prefix,
                       "write_frequency":write_frequency}

        self._prepare_calc(output_directory=output_directory,
                           calc_params=calc_params)
        
        f = open("output_file.txt","w")
        f.write("fake\n")
        f.close()

        self._complete_calc()

def test_Simulation(ens_test_data,variable_types):

    ens = ens_test_data["ens"]
    ddg_df = ens_test_data["ddg_df"]
    conditions = ens_test_data["conditions"]

    # Make sure implementation checks are in place. 
    with pytest.raises(NotImplementedError):
        sm = Simulation(ens=ens,
                        ddg_df=ddg_df,
                        conditions=conditions,
                        seed=None)


    with pytest.raises(NotImplementedError):
        sm = Simulation_no_run(ens=ens,
                               ddg_df=ddg_df,
                               conditions=conditions,
                               seed=None)
   
    with pytest.raises(NotImplementedError):
        sm = Simulation_no_calc_type(ens=ens,
                                     ddg_df=ddg_df,
                                     conditions=conditions,
                                     seed=None)

    # Now test variable checking

    # ens
    for v in variable_types["everything"]:
        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            sm = SimulationTester(ens=v,
                                  ddg_df=ddg_df,
                                  conditions=conditions,
                                  seed=None)
        
    # ddg_df
    for v in variable_types["everything"]:

        # skip strings -- leads to FileNotFound error because attempt to read
        # str as file
        if issubclass(type(v),str):
            continue

        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            sm = SimulationTester(ens=ens,
                                  ddg_df=v,
                                  conditions=conditions,
                                  seed=None)

    # conditions
    for v in variable_types["everything"]:

        # str will through FileNotFound, ignore
        if issubclass(type(v),str):
            continue

        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
    
            sm = SimulationTester(ens=ens,
                                  ddg_df=ddg_df,
                                  conditions=v,
                                  seed=None)

    # seed
    for v in variable_types["everything"]:
        
        # Skip coercable to float values
        try:
            int_v = int(v)
            if int_v >= 0:
                continue
        except:
            pass

        if v is None:
            continue

        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            sm = SimulationTester(ens=ens,
                                  ddg_df=ddg_df,
                                  conditions=conditions,
                                  seed=v)

    # Now test that things are being set reasonably well

    # ------------------------------------------------------------------
    # Check seed
    
    sm = SimulationTester(ens=ens,
                          ddg_df=ddg_df,
                          conditions=conditions,
                          seed=None)
    
    assert issubclass(type(sm._seed),int)
    assert sm._seed >= 0
    assert issubclass(type(sm._pcg64),np.random._pcg64.PCG64)

    sm = SimulationTester(ens=ens,
                          ddg_df=ddg_df,
                          conditions=conditions,
                          seed=5)

    assert issubclass(type(sm._seed),int)
    assert sm._seed == 5
    assert issubclass(type(sm._pcg64),np.random._pcg64.PCG64)
    assert np.isclose(0.8050029237453802,sm._rng.random())

    assert issubclass(type(sm._fc),Fitness)
    assert issubclass(type(sm._gc),Genotype)

def test_Simulation__prepare_calc(ens_test_data,tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    assert not os.path.exists("test_dir")

    ens = ens_test_data["ens"]
    ddg_df = ens_test_data["ddg_df"]
    conditions = ens_test_data["conditions"]

    sm = SimulationTester(ens=ens,
                          ddg_df=ddg_df,
                          conditions=conditions,
                          seed=5)

    sm._prepare_calc(output_directory="test_dir",
                     calc_params={"test":1})
    assert os.path.split(os.getcwd())[-1] == "test_dir"
    with open(os.path.join("input","simulation.json")) as f:
        sim_json = json.load(f)
    assert sim_json["calc_params"]["test"] == 1

    sm._complete_calc()
    assert os.getcwd() == tmpdir
    assert os.path.exists("test_dir")

    # Should not work now because directory exists
    with pytest.raises(FileExistsError):
        sm._prepare_calc(output_directory="test_dir",
                         calc_params={"test":1})


    os.chdir(current_dir)

def test_Simulation__write_calc_params(ens_test_data,
                                       tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    ens = ens_test_data["ens"]
    ddg_df = ens_test_data["ddg_df"]
    conditions = ens_test_data["conditions"]

    sm = SimulationTester(ens=ens,
                          ddg_df=ddg_df,
                          conditions=conditions,
                          seed=5)

    
    calc_params = {"population_size":100,
                   "mutation_rate":0.01,
                   "num_generations":100,
                   "write_prefix":"eee_sim",
                   "write_frequency":1000}
    
    sm._write_calc_params(calc_params=calc_params)

    assert os.path.exists(os.path.join("input","simulation.json"))
    assert os.path.exists(os.path.join("input","ddg.csv"))
    assert os.path.exists(os.path.join("input","ensemble.csv"))
    assert os.path.exists(os.path.join("input","conditions.csv"))

    with open(os.path.join("input","simulation.json")) as f:
        as_written = json.load(f)

    # Make sure run params written correctly. 
    for k in calc_params:
        assert as_written["calc_params"][k] == calc_params[k]
    
    # ens
    ens_written = read_ensemble(os.path.join("input",as_written["ens"])).to_dict()
    ens_dict = ens.to_dict()

    for k in ens_dict["ens"]:
        for n in ens_dict["ens"][k]:
            print(ens_dict["ens"][k])
            print(ens_written["ens"][k])
            assert ens_dict["ens"][k][n] == ens_written["ens"][k][n]

    assert ens_dict["gas_constant"] == as_written["gas_constant"]

    # conditions
    assert as_written["conditions"] == "conditions.csv"
    cond_written, _ = read_fitness_conditions(conditions=os.path.join("input",
                                                              "conditions.csv"),
                                      ens=sm.ens)
    
    keys_to_check =  ["select_on",
                      "select_on_folded",
                      "fitness_kwargs",
                      "temperature",
                      "fitness_fcn",
                      "X",
                      "Y"]
    for key in keys_to_check:
        assert np.array_equal(cond_written[key],sm.fc.condition_df[key])

    # ddg_df
    assert as_written["ddg_df"] == "ddg.csv"
    
    assert as_written["seed"] == sm._seed

    # Make sure we can read in dataframe
    df = pd.read_csv(os.path.join("input","ddg.csv"))

    os.chdir(current_dir)

def test_Simulation__complete_calc():
    # Tested within test__prepare_calc because these are paired functions. 
    assert True
    
def test_Simulation_system_params(ens_test_data):
        
    ens = ens_test_data["ens"]
    ddg_df = ens_test_data["ddg_df"]
    conditions = ens_test_data["conditions"]

    sm = SimulationTester(ens=ens,
                          ddg_df=ddg_df,
                          conditions=conditions,
                          seed=5)
    

    system_params = sm.system_params

    ens = sm._ens.to_dict()
    for k in system_params["ens"]:
        assert ens["ens"][k] == system_params["ens"][k]

    fc = sm._fc.to_dict()
    for k in fc:
        if hasattr(fc[k],"__iter__"):
            assert np.array_equal(fc[k],system_params["conditions"][k])
        else:
            assert fc[k] == system_params[k]

    gc = sm._gc.to_dict()
    assert gc["ddg_df"] is sm._gc._ddg_df

    assert system_params["seed"] == sm._seed


def test_Simulation_get_calc_description(ens_test_data):
    
    # print not a great test... mostly just make sure it runs without error and
    # produces a string. 

    ens = ens_test_data["ens"]
    ddg_df = ens_test_data["ddg_df"]
    conditions = ens_test_data["conditions"]

    sm = SimulationTester(ens=ens,
                          ddg_df=ddg_df,
                          conditions=conditions,
                          seed=5)
    
    # No kwargs passed in 
    value1 = sm.get_calc_description()
    assert issubclass(type(value1),str)

    # Kwargs passed in
    value2 = sm.get_calc_description(calc_kwargs={"mutation_rate":0.1})
    assert issubclass(type(value2),str)

    # kwargs should make longer
    assert len(value1) < len(value2)


def test_Simulation_fitness_from_energy(ens_test_data):
    

    ens = ens_test_data["ens"]
    ddg_df = ens_test_data["ddg_df"]
    conditions = ens_test_data["conditions"]

    sm = SimulationTester(ens=ens,
                          ddg_df=ddg_df,
                          conditions=conditions)
    
    f = sm.fitness_from_energy({"s1":0,"s2":0})
    assert f == np.prod(sm._gc._fitness_function(np.zeros(2)))

    f = sm.fitness_from_energy(np.zeros(2))
    assert f == np.prod(sm._gc._fitness_function(np.zeros(2)))

    
    mut_energy = np.array([1,0])
    f = sm.fitness_from_energy({"s1":1,"s2":0})
    assert f == np.prod(sm._gc._fitness_function(mut_energy))

    f = sm.fitness_from_energy(mut_energy)
    assert f == np.prod(sm._gc._fitness_function(mut_energy))

def test_Simulation_ens(ens_test_data):
    
    ens = ens_test_data["ens"]
    ddg_df = ens_test_data["ddg_df"]
    conditions = ens_test_data["conditions"]

    sm = SimulationTester(ens=ens,
                          ddg_df=ddg_df,
                          conditions=conditions)
    
    assert ens is sm.ens
    assert ens is sm._ens

def test_Simulation_fc(ens_test_data):
    
    ens = ens_test_data["ens"]
    ddg_df = ens_test_data["ddg_df"]
    conditions = ens_test_data["conditions"]

    sm = SimulationTester(ens=ens,
                          ddg_df=ddg_df,
                          conditions=conditions)
    
    assert sm.fc is sm._fc

def test_Simulation_gc(ens_test_data):
    
    ens = ens_test_data["ens"]
    ddg_df = ens_test_data["ddg_df"]
    conditions = ens_test_data["conditions"]

    sm = SimulationTester(ens=ens,
                          ddg_df=ddg_df,
                          conditions=conditions)
    
    assert sm.gc is sm._gc