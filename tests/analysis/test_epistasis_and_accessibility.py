
import pytest

from eee.analysis.epistasis_and_accessibility import _get_rings
from eee.analysis.epistasis_and_accessibility import epistasis_and_accessibility
from eee.core.data import AA_1TO3

import numpy as np
import pandas as pd

import os
import shutil
import json

def test__get_rings(tiny_sim_output,tmpdir):
    
    current_dir = os.getcwd()
    os.chdir(tmpdir)

    df, ring_one_muts, ring_two_muts, ring_two_genotypes = _get_rings(os.path.join(tiny_sim_output["dms"],
                                                                                   "eee_dms.csv"))
    assert len(df) == 400

    aa = [a for a in AA_1TO3.keys() if a != "L"]
    expected = []
    for a in aa:
        expected.append(f"L1{a}")
        expected.append(f"L2{a}")

    genotypes = []
    for a in aa:
        for b in aa:
            genotypes.append(f"L1{a}/L2{b}")

    assert set(expected) == ring_one_muts
    assert set(expected) == ring_two_muts
    assert set(genotypes) == ring_two_genotypes


    os.chdir(current_dir)

def test_epistasis_and_accessibility(tiny_sim_output,tmpdir):
    
    current_dir = os.getcwd()
    os.chdir(tmpdir)

    df = epistasis_and_accessibility(dms_dir=tiny_sim_output["dms"],
                                     accessible_dir=tiny_sim_output["acc"])
    
    assert issubclass(type(df),pd.DataFrame)
    assert len(df) == 7
    assert np.sum(df["ep_class"] == "mag") == 2
    assert np.sum(df["ep_class"] == "recip") == 5
    assert np.sum(df["stuck"]) == 3

    with pytest.raises(ValueError):
        df = epistasis_and_accessibility(dms_dir=tiny_sim_output["acc"],
                                         accessible_dir=tiny_sim_output["dms"])
        
    with pytest.raises(ValueError):
        df = epistasis_and_accessibility(dms_dir=tiny_sim_output["dms"],
                                         accessible_dir=tiny_sim_output["dms"])
        
    # ------------  diff ens

    # Copy in dms and acc
    shutil.copytree(tiny_sim_output["dms"],"dms")
    shutil.copytree(tiny_sim_output["acc"],"acc")

    # Should run
    df = epistasis_and_accessibility(dms_dir="dms",
                                     accessible_dir="acc")

    # modify ensemble so acc is different between two sims
    with open(os.path.join("acc","input","simulation.json")) as f:
        this_json = json.load(f)
    
    this_json["ens"]["hdna"]["dG0"] = 5

    with open(os.path.join("acc","input","simulation.json"),"w") as f:
        json.dump(this_json,f)
    
    with pytest.raises(ValueError):
        df = epistasis_and_accessibility(dms_dir="dms",
                                         accessible_dir="acc")
        
    # ------------  diff fitness

    # Copy in dms and acc
    shutil.rmtree("dms")
    shutil.rmtree("acc")
    shutil.copytree(tiny_sim_output["dms"],"dms")
    shutil.copytree(tiny_sim_output["acc"],"acc")

    # Should run
    df = epistasis_and_accessibility(dms_dir="dms",
                                     accessible_dir="acc")

    # modify fc so acc is different between two sims
    with open(os.path.join("acc","input","simulation.json")) as f:
        this_json = json.load(f)

    this_json["conditions"]["select_on_folded"] = [False,False]

    with open(os.path.join("acc","input","simulation.json"),"w") as f:
        json.dump(this_json,f)
    
    with pytest.raises(ValueError):
        df = epistasis_and_accessibility(dms_dir="dms",
                                         accessible_dir="acc")


    # ------------  acc max_depth parameter

    # Copy in dms and acc
    shutil.rmtree("dms")
    shutil.rmtree("acc")
    shutil.copytree(tiny_sim_output["dms"],"dms")
    shutil.copytree(tiny_sim_output["acc"],"acc")

    # Should run
    df = epistasis_and_accessibility(dms_dir="dms",
                                     accessible_dir="acc")
    
    
    with open(os.path.join("acc","input","simulation.json")) as f:
        this_json = json.load(f)
    
    this_json["calc_params"]["max_depth"] = 1

    with open(os.path.join("acc","input","simulation.json"),"w") as f:
        json.dump(this_json,f)
    
    with pytest.raises(ValueError):
        df = epistasis_and_accessibility(dms_dir="dms",
                                         accessible_dir="acc") 


    # ------------  dms max_depth parameter

    # Copy in dms and acc
    shutil.rmtree("dms")
    shutil.rmtree("acc")
    shutil.copytree(tiny_sim_output["dms"],"dms")
    shutil.copytree(tiny_sim_output["acc"],"acc")

    # Should run
    df = epistasis_and_accessibility(dms_dir="dms",
                                     accessible_dir="acc")
    
    
    with open(os.path.join("dms","input","simulation.json")) as f:
        this_json = json.load(f)
    
    this_json["calc_params"]["max_depth"] = 1

    with open(os.path.join("dms","input","simulation.json"),"w") as f:
        json.dump(this_json,f)
    
    with pytest.raises(ValueError):
        df = epistasis_and_accessibility(dms_dir="dms",
                                         accessible_dir="acc") 





    os.chdir(current_dir)