"""
Simulation class that finds all accessible pathways through sequence space. 
"""
from .simulation_base import Simulation
from eee.core.engine import pathfinder

from eee._private.check.standard import check_int
from eee._private.check.standard import check_bool
from eee._private.interface import run_cleanly

class AccessiblePaths(Simulation):

    calc_type = "accessible_paths"

    @run_cleanly
    def run(self,
            output_directory="eee_accessible",
            max_depth=1,
            allow_neutral=True,
            find_all_paths=True,
            output_file="eee_accessible.csv"):
        """
        Identify all accessible evolutionary paths starting from a wildtype
        protein.
        
        Parameters
        ----------
        output_directory : str, default="eee_sim"
            do simulation in this output directory
        max_depth : int, default=1
            explore paths up to this length. 1 corresponds to all accessible
            single mutants, 2 to all double mutants, 3 to all triple, etc. 
        output_file : str, default="eee_dms.csv"
            write results to the indicated csv file
        allow_neutral : bool, default=True
            allow mutations that have no effect on fitness
        find_all_paths : bool, default=True
            visit the same sequence from different starting points. This will find 
            all allowed paths to a given genotype. If False, ignore new paths that
            visit the same genotype. 
        output_file : str, default="genotypes.csv"
            return visited genotypes (with trajectory information) to this file
        """

        max_depth = check_int(value=max_depth,
                              variable_name="max_depth",
                              minimum_allowed=1)
        allow_neutral = check_bool(allow_neutral,
                                   variable_name="allow_neutral")
        find_all_paths = check_bool(value=find_all_paths,
                                    variable_name="find_all_paths")
        
        output_file = f"{output_file}"
    
        # Record the new keys
        calc_params = {}
        calc_params["max_depth"] = max_depth
        calc_params["allow_neutral"] = allow_neutral
        calc_params["find_all_paths"] = find_all_paths
        calc_params["output_file"] = output_file

        self._prepare_calc(output_directory=output_directory,
                           calc_params=calc_params)
        
        # Find all paths accessible from the wildtype genotype that involve 
        # single mutations and do not compromise fitness
        pathfinder(gc=self._gc,
                   max_depth=calc_params["max_depth"],
                   allow_neutral=calc_params["allow_neutral"],
                   find_all_paths=calc_params["find_all_paths"],
                   output_file=calc_params["output_file"],
                   return_output=False)
        
        self._complete_calc()


    