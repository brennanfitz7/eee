import pytest

from eee._private.check.eee.temperature import check_temperature

import numpy as np

def test_check_temperature(variable_types):
    
    allowed = variable_types["floats_or_coercable"][:]
    not_allowed = variable_types["not_floats_or_coercable"][:]
    for v in allowed:
        print(v,type(v),flush=True)

        if float(v) <= 0:
            not_allowed.append(v)
            continue

        out = check_temperature(temperature=v,num_conditions=1)
        assert len(out) == 1
        assert out[0] == float(v)

    for v in not_allowed:
        print(v,type(v),flush=True)

        with pytest.raises(ValueError):
            check_temperature(temperature=v,num_conditions=1)
    
    allowed = [[1,2],(1,2),np.array([1,2])]
    not_allowed = [["a","b"]]
    for v in allowed:
        out = check_temperature(temperature=v,num_conditions=2)
        assert len(out) == 2
        assert out[0] == float(v[0])
        assert out[1] == float(v[1])

        with pytest.raises(ValueError):
            out = check_temperature(temperature=v,num_conditions=3)