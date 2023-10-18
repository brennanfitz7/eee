import pytest

from eee.calcs.accessible_paths import AccessiblePaths
from eee.calcs import read_json

import os

def test_AccessiblePaths(ens_test_data):
    
    ens = ens_test_data["ens"]
    ddg_df = ens_test_data["ddg_df"]
    conditions = ens_test_data["conditions"]
    
    ap = AccessiblePaths(ens=ens,
                        ddg_df=ddg_df,
                        conditions=conditions,
                        seed=None)
    
    assert ap.calc_type == "accessible_paths"
    
def test_AccessiblePaths_run(ens_test_data,tmpdir):

    # This class directly wraps engine.pathfinder. Heavy testing there. 
    # Basically make sure it is working here. 

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    ens = ens_test_data["ens"]
    ddg_df = ens_test_data["ddg_df"]
    conditions = ens_test_data["conditions"]
    
    ap = AccessiblePaths(ens=ens,
                        ddg_df=ddg_df,
                        conditions=conditions,
                        seed=None)

    ap.run(output_directory="test",
           max_depth=2,
           allow_neutral=True,
           find_all_paths=True,
           output_file="yo.csv")
    
    assert os.path.exists(os.path.join("test","yo.csv"))
    assert os.path.exists(os.path.join("test","input","ddg.csv"))
    assert os.path.exists(os.path.join("test","input","ensemble.csv"))
    assert os.path.exists(os.path.join("test","input","conditions.csv"))
    assert os.path.exists(os.path.join("test","input","simulation.json"))

    os.chdir('test')
    _, kwargs = read_json(os.path.join("input",'simulation.json'))
    assert kwargs["max_depth"] == 2
    assert kwargs["allow_neutral"] == True
    assert kwargs["find_all_paths"] == True
    assert kwargs["output_file"] == "yo.csv"

    os.chdir("..")

 
    os.chdir(current_dir)