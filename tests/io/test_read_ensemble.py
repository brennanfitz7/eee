import pytest

from eee.core.data import GAS_CONSTANT
from eee.io.read_ensemble import _search_for_key
from eee.io.read_ensemble import _spreadsheet_to_ensemble
from eee.io.read_ensemble import _json_to_ensemble
from eee.io.read_ensemble import _file_to_ensemble
from eee.io.read_ensemble import read_ensemble

import numpy as np
import pandas as pd

import os
import shutil
import json

def test__search_for_key():

    some_dict = {"test":{"this":{"out":1},
                         "stupid":{1:5},
                         "rocket":{"ship":2}
                        },
                  "x":[1,2,3]}

    out = _search_for_key(some_dict,"test")
    assert np.array_equal(out,["test"])

    out = _search_for_key(some_dict,"this")
    assert np.array_equal(out,["test","this"])

    out = _search_for_key(some_dict,"ship")
    assert np.array_equal(out,["test","rocket","ship"])

    out = _search_for_key(some_dict,"x")
    assert np.array_equal(out,["x"])

    out = _search_for_key(some_dict,"not_a_key")
    assert len(out) == 0
    
def test__spreadsheet_to_ensemble(ensemble_inputs):
    
    df = ensemble_inputs["s1s2_folded.xlsx"]

    # Basic read 
    ens = _spreadsheet_to_ensemble(df=df)
    assert np.array_equal(ens.species,["s1","s2"])

    assert ens._species_dict["s1"]["dG0"] == 0
    assert ens._species_dict["s1"]["X"] == 1
    assert ens._species_dict["s1"]["Y"] == 0
    assert ens._species_dict["s1"]["folded"] == True
    assert ens._species_dict["s1"]["observable"] == True
        
    assert ens._species_dict["s2"]["dG0"] == 10
    assert ens._species_dict["s2"]["X"] == 0
    assert ens._species_dict["s2"]["Y"] == 1
    assert ens._species_dict["s2"]["folded"] == True
    assert ens._species_dict["s2"]["observable"] == False

    assert ens._gas_constant == GAS_CONSTANT

    # One change to entry 
    df = ensemble_inputs["s1s2_unfolded.xlsx"]

    ens = _spreadsheet_to_ensemble(df=df,
                                   gas_constant=1.0)
    assert np.array_equal(ens.species,["s1","s2"])

    assert ens._species_dict["s1"]["dG0"] == 0
    assert ens._species_dict["s1"]["X"] == 1
    assert ens._species_dict["s1"]["Y"] == 0
    assert ens._species_dict["s1"]["folded"] == True
    assert ens._species_dict["s1"]["observable"] == True
        
    assert ens._species_dict["s2"]["dG0"] == 10
    assert ens._species_dict["s2"]["X"] == 0
    assert ens._species_dict["s2"]["Y"] == 1
    assert ens._species_dict["s2"]["folded"] == False
    assert ens._species_dict["s2"]["observable"] == False   
    
    assert ens._gas_constant == 1.0

    # Drop some columns
    df = ensemble_inputs["s1s2_nostoich.xlsx"]

    ens = _spreadsheet_to_ensemble(df=df)
    assert np.array_equal(ens.species,["s1","s2"])

    assert ens._species_dict["s1"]["dG0"] == 0
    assert ens._species_dict["s1"]["folded"] == True
    assert ens._species_dict["s1"]["observable"] == True
    assert len(ens._species_dict["s1"]) == 3
        
    assert ens._species_dict["s2"]["dG0"] == 10
    assert ens._species_dict["s2"]["folded"] == False
    assert ens._species_dict["s2"]["observable"] == False   
    assert len(ens._species_dict["s2"]) == 3
    
    # Totally minimal
    df = pd.DataFrame({"name":["s1","s2","s3"]})

    ens = _spreadsheet_to_ensemble(df=df)
    assert ens._species_dict["s1"]["dG0"] == 0
    assert ens._species_dict["s1"]["folded"] == True
    assert ens._species_dict["s1"]["observable"] == False
    assert len(ens._species_dict["s1"]) == 3
        
    assert ens._species_dict["s2"]["dG0"] == 0
    assert ens._species_dict["s2"]["folded"] == True
    assert ens._species_dict["s2"]["observable"] == False  
    assert len(ens._species_dict["s2"]) == 3

    assert ens._species_dict["s3"]["dG0"] == 0
    assert ens._species_dict["s3"]["folded"] == True
    assert ens._species_dict["s3"]["observable"] == False  
    assert len(ens._species_dict["s3"]) == 3

    # Missing something required (name only required, actually)
    df = pd.DataFrame({"dG0":[1,3]})
    with pytest.raises(ValueError):
        ens = _spreadsheet_to_ensemble(df=df)

def test__json_to_ensemble(sim_json,ensemble_inputs,tmpdir):
    
    with open(sim_json["dms.json"]) as f:
        test_json = json.load(f)

    ens = _json_to_ensemble(test_json)

    assert ens._gas_constant == 0.001987

    assert ens._species_dict["hdna"]["dG0"] == 0
    assert ens._species_dict["hdna"]["folded"] == True
    assert ens._species_dict["hdna"]["observable"] == True
    assert len(ens._species_dict["hdna"]) == 3

    assert ens._species_dict["h"]["dG0"] == 5
    assert ens._species_dict["h"]["folded"] == True
    assert ens._species_dict["h"]["observable"] == False
    assert len(ens._species_dict["h"]) == 3

    assert ens._species_dict["l2e"]["dG0"] == 5
    assert ens._species_dict["l2e"]["iptg"] == 4
    assert ens._species_dict["l2e"]["folded"] == True
    assert ens._species_dict["l2e"]["observable"] == False
    assert len(ens._species_dict["l2e"]) == 4

    assert ens._species_dict["unfolded"]["dG0"] == 10
    assert ens._species_dict["unfolded"]["folded"] == False
    assert ens._species_dict["unfolded"]["observable"] == False
    assert len(ens._species_dict["unfolded"]) == 3

    with open(sim_json["lac.json"]) as f:
        test_json = json.load(f)

    ens = _json_to_ensemble(test_json)

    # Test gas constant reading

    assert ens._gas_constant == 0.008314

    assert ens._species_dict["hdna"]["dG0"] == 0
    assert ens._species_dict["hdna"]["folded"] == True
    assert ens._species_dict["hdna"]["observable"] == True
    assert len(ens._species_dict["hdna"]) == 3

    assert ens._species_dict["h"]["dG0"] == 5
    assert ens._species_dict["h"]["folded"] == True
    assert ens._species_dict["h"]["observable"] == False
    assert len(ens._species_dict["h"]) == 3

    assert ens._species_dict["l2e"]["dG0"] == 5
    assert ens._species_dict["l2e"]["iptg"] == 4
    assert ens._species_dict["l2e"]["folded"] == True
    assert ens._species_dict["l2e"]["observable"] == False
    assert len(ens._species_dict["l2e"]) == 4

    assert ens._species_dict["unfolded"]["dG0"] == 10
    assert ens._species_dict["unfolded"]["folded"] == False
    assert ens._species_dict["unfolded"]["observable"] == False
    assert len(ens._species_dict["unfolded"]) == 3

    # ens with no entries
    with open(ensemble_inputs["empty-ensemble.json"]) as f:
        test_json = json.load(f)
    ens = _json_to_ensemble(test_json)
    assert len(ens._species_dict) == 0
    assert ens._gas_constant == 0.008314

    # ensemble as top-level key
    with open(ensemble_inputs["top-level-ensemble.json"]) as f:
        test_json = json.load(f)
    ens = _json_to_ensemble(test_json)
    assert len(ens._species_dict) == 4
    assert ens._gas_constant == 0.008314

    # no gas constant defined
    with open(ensemble_inputs["no-gas-constant.json"]) as f:
        test_json = json.load(f)
    ens = _json_to_ensemble(test_json)
    assert len(ens._species_dict) == 4
    assert ens._gas_constant == GAS_CONSTANT

    # no ensemble key at all -- die.
    with open(ensemble_inputs["no-ensemble.json"]) as f:
        test_json = json.load(f)
    with pytest.raises(ValueError):
        ens = _json_to_ensemble(test_json) 

    with open(ensemble_inputs["ensemble-with-bad-key.json"]) as f:
        test_json = json.load(f)
    with pytest.raises(ValueError):
        ens = _json_to_ensemble(test_json) 

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    # Load ensemble json where ens has "spreadsheet" key and non-default
    # gas constant. 

    shutil.copy(ensemble_inputs["s1s2_folded.xlsx"],".")
    with open(ensemble_inputs["spreadsheet-ensemble.json"]) as f:
        test_json = json.load(f)
    ens = _json_to_ensemble(test_json)
    ens._gas_constant == 0.008314
    assert np.array_equal(ens.species,["s1","s2"])

    assert ens._species_dict["s1"]["dG0"] == 0
    assert ens._species_dict["s1"]["X"] == 1
    assert ens._species_dict["s1"]["Y"] == 0
    assert ens._species_dict["s1"]["folded"] == True
    assert ens._species_dict["s1"]["observable"] == True
        
    assert ens._species_dict["s2"]["dG0"] == 10
    assert ens._species_dict["s2"]["X"] == 0
    assert ens._species_dict["s2"]["Y"] == 1
    assert ens._species_dict["s2"]["folded"] == True
    assert ens._species_dict["s2"]["observable"] == False

    os.chdir(current_dir)

def test__file_to_ensemble(tmpdir,ensemble_inputs):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    with pytest.raises(FileNotFoundError):
        _file_to_ensemble("not_a_file")
    
    ens = _file_to_ensemble(ensemble_inputs["top-level-ensemble.json"])
    assert ens._gas_constant == 0.008314
    assert len(ens.species) == 4

    ens = _file_to_ensemble(ensemble_inputs["s1s2_folded.xlsx"])
    assert ens._gas_constant == GAS_CONSTANT
    assert len(ens.species) == 2

    os.chdir(current_dir)

def test_read_ensemble(tmpdir,ensemble_inputs,variable_types):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    # Send in something stupid
    with pytest.raises(FileNotFoundError):
        read_ensemble("not_a_file")

    # Send in json file
    ens = read_ensemble(ensemble_inputs["top-level-ensemble.json"])
    assert ens._gas_constant == 0.008314
    assert len(ens.species) == 4

    # Send in excel file
    ens = read_ensemble(ensemble_inputs["s1s2_folded.xlsx"])
    assert ens._gas_constant == GAS_CONSTANT
    assert len(ens.species) == 2
    
    # Send in pandas 
    df = pd.read_excel(ensemble_inputs["s1s2_folded.xlsx"])
    ens = read_ensemble(df)
    assert ens._gas_constant == GAS_CONSTANT
    assert len(ens.species) == 2

    # Send in json 
    with open(ensemble_inputs["top-level-ensemble.json"],'r') as f:
        json_input = json.load(f)
    
    ens = read_ensemble(json_input)
    assert ens._gas_constant == 0.008314
    assert len(ens.species) == 4 

    # General error checking
    for v in variable_types["everything"]:
        
        expected_err = ValueError
        if issubclass(type(v),dict):
            continue
        if issubclass(type(v),str):
            expected_err = FileNotFoundError
        
        if issubclass(type(v),pd.DataFrame):
            continue


        print(v,type(v))
        with pytest.raises(expected_err):
            read_ensemble(v)
        
    os.chdir(current_dir)