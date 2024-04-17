"""
Code for running evolutionary simulations. Main user-facing simulation code. 
All Simulation subclasses imported here will be available for users to access
in simulation.json files using the `calc_type` key.
"""

from .simulation_base import Simulation as _SimulationBaseClass

from .wf_sim import WrightFisherSimulation
from .dms import DeepMutationalScan
from .accessible_paths import AccessiblePaths
from .wf_tree_sim import WrightFisherTreeSimulation

def _get_available():

    calc_available = {}

    possible = dict(globals())
    for p in possible:
        this_poss = possible[p]
        
        # Get any subclass of Simulation base class and key to calc_type in 
        # calc_available
        try:
            if issubclass(this_poss,_SimulationBaseClass):
                if this_poss.calc_type is not None:
                    calc_available[this_poss.calc_type] = this_poss
        except TypeError:
            continue
    
    return calc_available
                
# Register calculations available
CALC_AVAILABLE = _get_available()

from .read_json import read_json
