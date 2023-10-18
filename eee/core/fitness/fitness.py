"""
Class for calculating fitness from an ensemble during an evolutionary simulation. 
"""

from .read_fitness_conditions import read_fitness_conditions
from .map_fitness_fcn import map_fitness_fcn

from eee._private.check.ensemble import check_ensemble

import numpy as np

import copy

class Fitness:
    """
    Class used to calculate fitness of genotypes in evolutionary simulations.
    It holds fixed aspects of the fitness (ensemble, ligand chemical potentials,
    fitness functions, and temperature), allowing us to quickly calculate
    fitness given only the energy of a particular genotype. 

    Notes
    -----
    This class stores a private *copy* of ens. This is so the ensemble z-matrix 
    stays identical for all calculations, even if the user uses the initial 
    ensemble object for a calculation with a different ligand_dict after initialization
    of the FitnessContainer object. (That calculation triggers creation of new
    z-matrix, which would potentially change output observables). 
    """
    def __init__(self,
                 ens,
                 conditions,
                 default_fitness_kwargs=None,
                 default_select_on="fx_obs",
                 default_select_on_folded=True,
                 default_temperature=298.15):
        """
        Parameters
        ----------
        ens : eee.core.Ensemble 
            initialized instance of an Ensemble class
        conditions : pandas.DataFrame or similar
            Conditions at which to do the fitness calculation. Columns are 
            parameters, rows are conditions under which to do the calculation. 
            The `fitness_fcn` column is required. This indicates which fitness
            function to apply at the particular condition. Options (at this 
            writing) are "on", "off", "neutral", "on_above", and "on_below." 
            Other columns are: 
            
             + `fitness_kwargs`: keywords to pass to `fitness_fcn` (for example,
               `{"threshold":0.5}` for `on_above` and `on_below`). 
             + `select_on`: "fx_obs" or "dG_obs". All rows must have the same 
               value. 
             + `select_on_folded`: (True or False).
             + `temperature`: (temperature in K).

            All other columns are interpreted as ligand concentrations. The 
            column names much match ligands defined in `ens`. 
        default_fitness_kwargs : dict, optional
            if fitness_kwargs are not specified in conditions, assign this value
        default_select_on : str, default="fx_obs"
            if select_on is not specified in conditions, assign this value
        default_select_on_folded : bool, default=True
            if select_on_folded is not specified in conditions, assign this value
        default_temperature : float, default=298.15
            if temperature is not specified in conditions, assign this value
            
        Notes
        -----
        The conditions keyword argument is quite flexible. It can take a fully
        filled out dataframe, the filename of a dataframe, a dictionary with 
        columns as keys and list-like values holding entries, or a list of 
        dictionaries with values for each condition. 
        """
        
        # Load and validate the ensemble
        self._ens = check_ensemble(ens,check_obs=True)

        # Load the conditions
        out = read_fitness_conditions(conditions=conditions,
                              ens=self._ens,
                              default_fitness_kwargs=default_fitness_kwargs,
                              default_select_on=default_select_on,
                              default_select_on_folded=default_select_on_folded,
                              default_temperature=default_temperature) 
                              
        self._condition_df = out[0]
        self._ligand_dict = out[1]
        
        # Convert parts of condition_df to numpy arrays for ease of access
        self._select_on_folded = np.array(self._condition_df["select_on_folded"])
        self._fitness_kwargs = np.array(self._condition_df["fitness_kwargs"])
        self._temperature = np.array(self._condition_df["temperature"])
            
        # Set up ensemble to do calculation
        self._private_ens = copy.deepcopy(self._ens)
        self._private_ens.read_ligand_dict(ligand_dict=self._ligand_dict)
        self._select_on = self._condition_df["select_on"].iloc[0]
        self._obs_function = self._private_ens.get_observable_function(self._select_on)

        # Set up fitness calculations
        fitness_fcns = []
        for ff in self._condition_df["fitness_fcn"]:
            fitness_fcns.append(map_fitness_fcn(value=ff,return_as="function"))
        self._fitness_fcns = np.array(fitness_fcns)
                                
        self._num_conditions = len(self._fitness_fcns)
        self._F_array = np.zeros(self._num_conditions,dtype=float)


    def fitness(self,mut_energy_array):
        """
        Calculate the fitness of a genotype with total mutational energies 
        given by mut_energy_array. Fitness is defined as the product of the
        fitness  in each of the conditions specified in ligand_dict. 

        mut_energy_array : numpy.ndarray
            array holding the effects of mutations on energy. values should be
            in the order of ens.species 
        """
        
        values, fx_folded = self._obs_function(mut_energy_array=mut_energy_array,
                                               temperature=self._temperature)

        for i in range(self._num_conditions):
            self._F_array[i] = self._fitness_fcns[i](values[i],
                                                     **self._fitness_kwargs[i])
            if self._select_on_folded[i]:
               self._F_array[i] *= fx_folded[i]
 
        return self._F_array
    
    def to_dict(self):
        """
        Return a json-able dictionary describing the fitness parameters.
        """

        out = {}
        for k in self._condition_df.keys():
            out[k] = np.array(self._condition_df[k])
        return out

    @property
    def ens(self):
        return self._ens
        
    @property
    def ligand_dict(self):
        return self._ligand_dict
    
    @property
    def select_on(self):
        return self._select_on

    @property
    def select_on_folded(self):
        return self._select_on_folded

    @property
    def fitness_kwargs(self):
        return self._fitness_kwargs    
    
    @property
    def temperature(self):
        return self._temperature
    
    @property
    def condition_df(self):
        """
        Conditions as a pandas dataframe. 
        """

        return self._condition_df.copy()

        