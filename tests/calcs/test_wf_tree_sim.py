import pytest

from eee.calcs import WrightFisherTreeSimulation
from eee.calcs import read_json

import ete3

import os

def test_WrightFisherTreeSimulation(ens_test_data):
    
    ens = ens_test_data["ens"]
    ddg_df = ens_test_data["ddg_df"]
    conditions = ens_test_data["conditions"]

    wf = WrightFisherTreeSimulation(ens=ens,
                                    ddg_df=ddg_df,
                                    conditions=conditions, 
                                    seed=None)
    
    assert wf.calc_type == "wf_tree_sim"
    
def test_WrightFisherTreeSimulation_run(ens_test_data,
                                        newick_files,
                                        tmpdir):

    # Make sure wrapper runs. 

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    ens = ens_test_data["ens"]
    ddg_df = ens_test_data["ddg_df"]
    conditions = ens_test_data["conditions"]
    
    wf = WrightFisherTreeSimulation(ens=ens,
                                    ddg_df=ddg_df,
                                    conditions=conditions,
                                    seed=None)

    wf.run(output_directory="test",
           tree=ete3.Tree(newick_files["simple.newick"]),
           population_size=100,
           mutation_rate=0.1,
           num_generations=1000,
           burn_in_generations=10,
           write_prefix="eee_sim")
    
    expected_files = ["eee_sim.newick",
                      "eee_sim_genotypes.csv",
                      "eee_sim_burn-in-anc00.pickle",
                      "eee_sim_anc00-anc01.pickle",
                      "eee_sim_anc01-A.pickle",
                      "eee_sim_anc01-B.pickle",
                      "eee_sim_anc00-anc02.pickle",
                      "eee_sim_anc02-C.pickle",
                      "eee_sim_anc02-D.pickle"]
    
    for f in expected_files:
        assert os.path.exists(os.path.join("test",f))
 
    assert os.path.exists(os.path.join("test","input","ddg.csv"))
    assert os.path.exists(os.path.join("test","input","ensemble.csv"))
    assert os.path.exists(os.path.join("test","input","conditions.csv"))
    assert os.path.exists(os.path.join("test","input","simulation.json"))

    # This test makes sure the variables are being set properly and then 
    # recorded into the json. 
    os.chdir('test')
    _, kwargs = read_json(os.path.join("input",'simulation.json'))
    assert issubclass(type(kwargs["tree"]),ete3.TreeNode)
    assert kwargs["population_size"] == 100
    assert kwargs["mutation_rate"] == 0.1
    assert kwargs["num_generations"] == 1000
    assert kwargs["burn_in_generations"] == 10

    os.chdir("..")

    os.chdir(current_dir)

