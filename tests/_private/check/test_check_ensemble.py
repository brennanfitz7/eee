
from eee._private.check.ensemble import check_ensemble
from eee.core.ensemble import Ensemble

import pytest

def test_check_ensemble(variable_types):
    
    for v in variable_types:
        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            check_ensemble(v)

    ens = Ensemble(gas_constant=1)
    check_ensemble(ens)
    check_ensemble(ens,check_obs=False)

    with pytest.raises(ValueError):
        check_ensemble(ens,check_obs=True)

    ens.add_species(name="test1",
                    observable=True,
                    X=1)

    with pytest.raises(ValueError):
        check_ensemble(ens,check_obs=True)

    ens.add_species(name="test2",
                    observable=False,
                    Y=1)
    
    check_ensemble(ens,check_obs=True)