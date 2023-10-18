import pytest

from eee._private.check.eee.ligand_dict import check_ligand_dict

import pandas as pd
import numpy as np


def test_check_ligand_dict(variable_types):
        
    for v in variable_types["not_dict"]:
        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            check_ligand_dict(v)

    out, length = check_ligand_dict({})
    assert issubclass(type(out),dict)
    assert len(out) == 0
    assert length == 1

    bad_types = variable_types["types"][:]
    bad_types.extend([pd.DataFrame({"a":[1,2,3]}),
                      {"a":[1,2,3]}])
    for v in bad_types:
        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            check_ligand_dict({"X":v})

    bad_values = [[],np.array([])]
    for v in bad_values:
        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            check_ligand_dict({"X":v})

    bad_values = ["should_be_a_number"]
    for v in bad_values:
        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            check_ligand_dict({"X":v})

    good_values = variable_types["floats_or_coercable"][:]
    good_values.extend([[1,2,3],
                        [1.0,2.0,3.0],
                        [1,2.0,3],
                        np.array([1,2,3],dtype=int),
                        np.array([1,2,3],dtype=float)])
    for v in good_values:
        print(v,type(v),flush=True)
        check_ligand_dict({"X":v})
    
    good_values = [{"X":1,"Y":1},
                   {"X":"1","Y":"1"},
                   {"X":1,"Y":[1]},
                   {"X":[1],"Y":[1]}]
    for v in good_values:
        print(v,type(v),flush=True)
        result, length = check_ligand_dict(v)
        assert length == 1
        assert np.array_equal(result["X"],[1])
        assert np.array_equal(result["Y"],[1])

    good_values = [{"X":1,"Y":[1,2,3]},
                   {"X":[1],"Y":[1,2,3]},
                   {"X":[1,1,1],"Y":[1,2,3]}]
    for v in good_values:
        print(v,type(v),flush=True)
        result, length = check_ligand_dict(v)
        assert length == 3
        assert np.array_equal(result["X"],[1,1,1])
        assert np.array_equal(result["Y"],[1,2,3])

    good_values = [{"X":1,"Y":[1,2,3],"Z":[2]},
                   {"X":[1],"Y":[1,2,3],"Z":[2,2,2]},
                   {"X":[1,1,1],"Y":[1,2,3],"Z":[2,2,2]}]
    for v in good_values:
        print(v,type(v),flush=True)
        result, length = check_ligand_dict(v)
        assert length == 3
        assert np.array_equal(result["X"],[1,1,1])
        assert np.array_equal(result["Y"],[1,2,3])
        assert np.array_equal(result["Z"],[2,2,2])

    bad_values = [{"X":[1,2],"Y":[1,2,3]},
                  {"X":[1,2],"Y":[1,2,3],"Z":2}]
    for v in bad_values:
        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            check_ligand_dict(v)