"""
Code to run a Wright-Fisher simulation on an ensemble following an evolutionary
tree. 
"""
from eee.simulation.core.simulation import Simulation

from eee.simulation.core.engine import follow_tree

from eee._private.check.eee import check_num_generations
from eee._private.check.eee import check_mutation_rate
from eee._private.check.eee import check_population_size
from eee._private.check.eee import check_burn_in_generations
from eee._private.interface import run_cleanly


class WrightFisherTreeSimulation(Simulation):

    calc_type = "wf_tree_sim"

    @run_cleanly
    def run(self,
            newick,
            output_directory="eee_wf-tree-sim",
            population_size=1000,
            mutation_rate=0.01,
            num_generations=100,
            burn_in_generations=100,
            write_prefix="eee_wf-tree-sim"):
        """
        Run a simulation and save files to an output directory.
        
        Parameters
        ----------
        newick : str or ete3.Tree
            newick formatted tree with branch lengths and tip labels
        output_directory : str, default="eee_sim"
            do simulation in this output directory
        population_size : int, default=1000
            population size for the simulation. Should be > 0.
        mutation_rate : float, default=0.01
            mutation rate for the simulation. Should be > 0.
        num_generations : int, default=100
            this specifies the maximum number of generations allowed on each 
            branch. The simulation will run along each branch until either the
            correct number of mutations have accumulated (for branch length) 
            or the simulation hits num_generations.
        burn_in_generations : int, default=10
            run a Wright-Fisher simulation burn_in_generations long to generate an
            ancestral population. Must be >= 0. 
        write_prefix : str, default="eee_sim"
            write output files during the run with this prefix. 
        """

        population_size = check_population_size(population_size)
        mutation_rate = check_mutation_rate(mutation_rate) 
        num_generations = check_num_generations(num_generations)
        burn_in_generations = check_burn_in_generations(burn_in_generations)
        write_prefix = f"{write_prefix}"
    
        # Record the new keys
        calc_params = {}
        calc_params["newick"] = f"{write_prefix}.newick"
        calc_params["population_size"] = population_size
        calc_params["mutation_rate"] = mutation_rate
        calc_params["num_generations"] = num_generations
        calc_params["write_prefix"] = write_prefix
        calc_params["burn_in_generations"] = burn_in_generations

        self._prepare_calc(output_directory=output_directory,
                           calc_params=calc_params)
        
        # Run and return a Wright Fisher simulation.
        self._gc, _ =  follow_tree(gc=self._gc,
                                   newick=newick, 
                                   population=calc_params["population_size"],
                                   mutation_rate=calc_params["mutation_rate"],
                                   num_generations=calc_params["num_generations"],
                                   burn_in_generations=calc_params["burn_in_generations"],
                                   write_prefix=calc_params["write_prefix"],
                                   rng=self._rng)
        
        self._complete_calc()


    