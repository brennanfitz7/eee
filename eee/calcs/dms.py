"""
Run a deep mutational scan on a protein. 
"""
from .simulation_base import Simulation
from eee.core.engine import exhaustive

from eee._private.check.standard import check_int
from eee._private.interface import run_cleanly


class DeepMutationalScan(Simulation):

    calc_type = "dms"

    @run_cleanly
    def run(self,
            output_directory="eee_dms",
            max_depth=1,
            output_file="eee_dms.csv"):
        """
        Run a deep mutational scan up to max_depth mutations away from wildtype. 
        
        Parameters
        ----------
        output_directory : str, default="eee_sim"
            do simulation in this output directory
        max_depth : int, default=1
            max_depth for the deep-mutational scan. 1 corresponds to all single mutants,
            2 to all double mutants, 3 to all triple, etc. WARNING: The space gets
            very large as the number of sites and number of possible mutations 
            increase. 
        output_file : str, default="eee_dms.csv"
            write results to the indicated csv file
        """

        max_depth = check_int(value=max_depth,
                              variable_name="max_depth",
                              minimum_allowed=0)
        output_file = f"{output_file}"
    
        # Record the new keys
        calc_params = {}
        calc_params["max_depth"] = max_depth
        calc_params["output_file"] = output_file

        self._prepare_calc(output_directory=output_directory,
                           calc_params=calc_params)
        
        # Run and return a Wright Fisher simulation.
        exhaustive(gc=self._gc,
                   max_depth=calc_params["max_depth"],
                   output_file=calc_params["output_file"],
                   return_output=False)
        
        self._complete_calc()


    