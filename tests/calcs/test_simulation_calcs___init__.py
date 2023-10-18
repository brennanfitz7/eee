
from eee.calcs.simulation_base import Simulation
from eee.calcs import CALC_AVAILABLE

def test__get_calc_available():

    # Run on init, populating CALC_AVAILABLE

    assert len(CALC_AVAILABLE) > 0
    assert issubclass(type(CALC_AVAILABLE),dict)
    
    for k in CALC_AVAILABLE:
        assert issubclass(type(k),str)
        assert issubclass(CALC_AVAILABLE[k],Simulation)