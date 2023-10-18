import pytest

from eee._private.utility import prep_for_json

import numpy as np
import json

def test_prep_for_json(variable_types):

    # Make sure type check works
    for v in variable_types["everything"]:
        if issubclass(type(v),dict):
            continue
        print(v,type(v))
        with pytest.raises(ValueError):
            prep_for_json(v)

    # Test work on copy flag
    some_dict = {}
    out_dict = prep_for_json(some_dict=some_dict,
                             work_on_copy=True)
    assert out_dict is not some_dict

    out_dict = prep_for_json(some_dict=some_dict,
                             work_on_copy=False)
    assert out_dict is some_dict

    # Single int
    tests = [("int",np.int16,int,5),
             ("float",np.float16,float,5),
             ("bool",np.bool_,bool,True)]
    
    for t in tests:

        print(f"testing single {t[0]}")
        np_type = t[1]
        python_type = t[2]
        value = t[3]

        test_dict = {"key":np_type(value)}
        assert issubclass(type(test_dict["key"]),np_type)
        with pytest.raises(TypeError):
            json.dumps(test_dict)
        

        out_dict = prep_for_json(some_dict=test_dict)
        assert not issubclass(type(out_dict["key"]),np_type)
        assert issubclass(type(out_dict["key"]),python_type)
        assert out_dict["key"] == value
        json.dumps(out_dict)


    for t in tests:
        
        print(f"testing array of {t[0]}")
        np_type = t[1]
        python_type = t[2]
        value = t[3]

        array = np.ones(5,dtype=np_type)*value
        test_dict["key"] = array
        assert issubclass(type(test_dict["key"]),np.ndarray)
        assert issubclass(type(test_dict["key"][0]),np_type)
        test_dict = {"key":array}
        with pytest.raises(TypeError):
            json.dumps(test_dict)

        out_dict = prep_for_json(some_dict=test_dict)
        assert not issubclass(type(out_dict["key"]),np.ndarray)
        assert not issubclass(type(out_dict["key"][0]),np_type)
        assert issubclass(type(out_dict["key"]),list)
        assert issubclass(type(out_dict["key"][0]),python_type)
        assert out_dict["key"][0] == value
        json.dumps(out_dict)

    # Send in dictionary with all kinds of wackiness and just test for 
    # writablity before and after
    array = np.ones(5,dtype=np.int16)*5
    nested_dict = { "l0.0":array.copy(),
                    "l0.1":{
                        "l1.0":array.copy(),
                        "l1.1":{
                            "l2.0":array.copy(),
                            "l2.1":np.int16(1),
                            "l2.2":np.float16(1),
                            "l2.3":np.bool_(True)},
                        "l1.2":np.int16(1),
                        "l1.3":np.float16(1),
                        "l1.4":np.bool_(True)},
                    "l0.2":np.int16(1),
                    "l0.3":np.float16(1),
                    "l0.4":np.bool_(True)}
    with pytest.raises(TypeError):
        json.dumps(nested_dict)
    
    nested_dict_out = prep_for_json(nested_dict)
    json.dumps(nested_dict_out)
    
    # Test mixture of bad and okay datatypes
    nested_dict = {"l0.0":array.copy(),
                   "l0.1":"test",
                   "l0.2":[1,2,3],
                   "l0.3":5,
                   "l0.4":np.bool_(False),
                   "l0.5":["test",1],
                   "l0.6":None}
    with pytest.raises(TypeError):
        json.dumps(nested_dict)
    nested_dict_out = prep_for_json(nested_dict)
    json.dumps(nested_dict_out)
    
    

