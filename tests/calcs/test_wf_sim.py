import pytest

from eee.calcs import WrightFisherSimulation
from eee.calcs import read_json

import pandas as pd

import os

def test_WrightFisherSimulation(ens_test_data):
    
    ens = ens_test_data["ens"]
    ddg_df = ens_test_data["ddg_df"]
    conditions = ens_test_data["conditions"]

    wf = WrightFisherSimulation(ens=ens,
                                ddg_df=ddg_df,
                                conditions=conditions,
                                seed=None)
    
    assert wf.calc_type == "wf_sim"
    
def test_WrightFisherSimulation_run(ens_test_data,tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    ens = ens_test_data["ens"]
    ddg_df = ens_test_data["ddg_df"]
    conditions = ens_test_data["conditions"]

    wf = WrightFisherSimulation(ens=ens,
                                ddg_df=ddg_df,
                                conditions=conditions,
                                seed=None)

    wf.run(output_directory="test",
           population_size=100,
           mutation_rate=0.01,
           num_generations=10000,
           num_mutations=1,
           write_prefix="eee_sim",
           write_frequency=1000,
           verbose=False)
    
    assert os.path.exists(os.path.join("test","input","ddg.csv"))
    assert os.path.exists(os.path.join("test","input","ensemble.csv"))
    assert os.path.exists(os.path.join("test","input","conditions.csv"))
    assert os.path.exists(os.path.join("test","input","simulation.json"))
    assert os.path.exists(os.path.join("test","eee_sim_genotypes.csv"))
    assert os.path.exists(os.path.join("test","eee_sim_generations_00.pickle"))
 
    os.chdir('test')
    
    _, kwargs = read_json(os.path.join("input",'simulation.json'))
    assert kwargs["population_size"] == 100
    assert kwargs["mutation_rate"] == 0.01
    assert kwargs["num_generations"] == 10000
    assert kwargs["num_mutations"] == 1
    assert kwargs["write_prefix"] == "eee_sim"
    assert kwargs["write_frequency"] == 1000
    assert kwargs["verbose"] == False

    os.chdir("..")

    os.chdir(current_dir)