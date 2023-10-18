import pytest

from eee.analysis.num_genotypes import get_num_genotypes

import numpy as np

def test_get_num_genotypes(variable_types):

    ddg_dict = {1:["A","B","C"],
                2:["E","F","G"]}

    shells = get_num_genotypes(ddg_dict=ddg_dict,
                               max_depth=0)
    assert np.array_equal(shells,[1])

    shells = get_num_genotypes(ddg_dict=ddg_dict,
                               max_depth=1)
    assert np.array_equal(shells,[1,6])

    shells = get_num_genotypes(ddg_dict=ddg_dict,
                               max_depth=2)
    assert np.array_equal(shells,[1,6,9])

    shells = get_num_genotypes(ddg_dict=ddg_dict,
                               max_depth=3)
    assert np.array_equal(shells,[1,6,9])


    ddg_dict = {1:["A","B","C"],
                2:["E","F","G"]}

    shells = get_num_genotypes(ddg_dict=ddg_dict,
                               max_depth=0)
    assert np.array_equal(shells,[1])

    shells = get_num_genotypes(ddg_dict=ddg_dict,
                               max_depth=1)
    assert np.array_equal(shells,[1,6])

    shells = get_num_genotypes(ddg_dict=ddg_dict,
                               max_depth=2)
    assert np.array_equal(shells,[1,6,9])

    shells = get_num_genotypes(ddg_dict=ddg_dict,
                               max_depth=3)
    assert np.array_equal(shells,[1,6,9])


    ddg_dict = {1:["A","B","C"],
                2:["E","F","G"],
                3:["H"]}
    
    shells = get_num_genotypes(ddg_dict=ddg_dict,
                               max_depth=0)
    assert np.array_equal(shells,[1])

    shells = get_num_genotypes(ddg_dict=ddg_dict,
                               max_depth=1)
    assert np.array_equal(shells,[1,7])

    shells = get_num_genotypes(ddg_dict=ddg_dict,
                               max_depth=2)
    assert np.array_equal(shells,[1,7,15])

    shells = get_num_genotypes(ddg_dict=ddg_dict,
                               max_depth=3)
    assert np.array_equal(shells,[1,7,15,9])

    # Send in bad stuff
    # ddg_dict
    for v in variable_types["everything"]:
        if issubclass(type(v),dict):
            continue

        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            get_num_genotypes(ddg_dict=v,
                              max_depth=2)

    # max_depth
    for v in variable_types["not_ints_or_coercable"]:
        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            get_num_genotypes(ddg_dict=ddg_dict,
                              max_depth=v)

    with pytest.raises(ValueError):
        get_num_genotypes(ddg_dict=ddg_dict,
                            max_depth=-1)
        
    
    
