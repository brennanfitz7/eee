"""
Function for calculating fitness from an ensemble.
"""


from eee.core.fitness import Fitness

def ensemble_fitness(ens,
                     conditions,
                     mut_energy=None):
    """
    Calculate fitness from the ensemble given mutations (in mut_energy), 
    chemical potentials (ligand_dict), and fitness function(s). 
    
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

    Returns
    -------
    df : pandas.DataFrame
        pandas dataframe holding information about the ensemble versus ligand
        concentration, including fitness. 
    """

    fc = Fitness(ens=ens,
                 conditions=conditions)
        
    df = ens.get_obs(mut_energy=mut_energy,
                     ligand_dict=fc.ligand_dict,
                     temperature=fc.temperature)

    mut_energy_array = ens.mut_dict_to_array(mut_energy=mut_energy)
    fitness_values = fc.fitness(mut_energy_array=mut_energy_array)
    df["fitness"] = fitness_values
    
    return df
