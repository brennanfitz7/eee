"""
Read eee conditions dataframe.
"""

from .map_fitness_fcn import map_fitness_fcn

from eee._private.check.ensemble import check_ensemble
from eee._private.check.dataframe import check_dataframe
from eee._private.check.standard import check_bool
from eee._private.check.standard import check_float

import json

def read_fitness_conditions(conditions,
                            ens,
                            default_fitness_kwargs=None,
                            default_select_on="fx_obs",
                            default_select_on_folded=True,
                            default_temperature=298.15):
    """
    Read a set of conditions for a simulation or fitness calculation. 

    Parameters
    ----------
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
    ens : eee.core.Ensemble 
        initialized instance of an Ensemble class
    default_fitness_kwargs : dict, optional
        if fitness_kwargs are not specified in conditions, assign this value
    default_select_on : str, default="fx_obs"
        if select_on is not specified in conditions, assign this value
    default_select_on_folded : bool, default=True
        if select_on_folded is not specified in conditions, assign this value
    default_temperature : float, default=298.15
        if temperature is not specified in conditions, assign this value

    Returns
    -------
    condition_df : pandas.DataFrame
        conditions in a standardized, validated pandas DataFrame
    ligand_dict : dict
        dictionary of ligands and their concentrations in a form that can be
        passed directly into an Ensemble object
    """

    # Read spreadsheet or list of conditions into a dataframe
    df = check_dataframe(conditions,
                         variable_name="conditions")
    df_columns = set(df.columns)

    # Set up defaults for values that are undefined in the conditions
    if default_fitness_kwargs is None:
        default_fitness_kwargs = {}

    defaults = {"fitness_kwargs":default_fitness_kwargs,
                "select_on":default_select_on,
                "select_on_folded":default_select_on_folded,
                "temperature":default_temperature}

    # Absolutely required columns
    required = set(["fitness_fcn"])

    # Reserved columns (built from defaults and required)
    reserved = list(defaults.keys())
    reserved.extend(list(required))
    reserved = set(reserved)
    
    # Make sure the dataframe has all required columns
    missing_required = list(required - df_columns)
    if len(missing_required) > 0:
        err = "\ndataframe is missing required columns:\n"
        for m in missing_required:
            err += f"    {m}\n"
        err += "\n"
        raise ValueError(err)
    
    # Get ligands that are defined in the ensemble
    ens = check_ensemble(ens,check_obs=False)
    ligands = set(ens.ligands)

    # Make sure the ensemble does not have any ligands with reserved names
    in_reserved = list(ligands.intersection(reserved))
    if len(in_reserved) > 0:
        err = "\nensemble must not have ligands with reserved names:\n"
        for r in in_reserved:
            err += f"    {r}\n"
        err += "\n"
        raise ValueError(err)
    
    # Make sure all non reserved columns correspond to ligands in the ensemble
    ligand_columns = (df_columns - reserved)
    extra_ligands = list(ligand_columns - ligands)
    if len(extra_ligands) > 0:
        err = "\ncolumns in the dataframe could not be interpreted as ligands.\n"
        for lig in extra_ligands:
            err += f"    {lig}\n"
        err += "\n"
        raise ValueError(err)

    # Get default values for columns not defined already
    for c in reserved:
        if c not in df_columns:
            df[c] = [defaults[c] for _ in range(len(df))]
    
    # Set ligands not in the conditions to 0.0
    for lig in ligands:
        if lig not in df_columns:
            df[lig] = 0.0

    # Parse and validate fitness functions
    fitness_fcn_out = []
    fitness_fcn_fcn = []
    for idx in df.index:
        v = df.loc[idx,"fitness_fcn"]
        f = map_fitness_fcn(v,return_as="function")
        v = map_fitness_fcn(f,return_as="string")
        fitness_fcn_out.append(v)
        fitness_fcn_fcn.append(f)
    df["fitness_fcn"] = fitness_fcn_out

    # fitness_kwargs
    kwargs_out = []
    for i, idx in enumerate(df.index):
        v = df.loc[idx,"fitness_kwargs"]

        # Convert into a dict if necessary
        if issubclass(type(v),str):
            v = json.loads(v)
        elif issubclass(type(v),dict):
            pass
        else:
            err = "\ncould not parse fitness_kwargs '{v}'\n\n"
            raise ValueError(err) 
        
        # Make sure we can use these fitness kwargs with the specified 
        # fitness function. 
        try:
            fitness_fcn_fcn[i](1.0,**v)
        except Exception as e:
            err = f"\nfitness_fcn {df.loc[idx,'fitness_fcn']} cannot take\n"
            err += f"fitness_kwargs '{df.loc[idx,'fitness_kwargs']}'\n\n"
            raise ValueError(err) from e

        kwargs_out.append(v)

    df["fitness_kwargs"] = kwargs_out
    
    # select_on
    select_on_out = []
    for idx in df.index:
        v = df.loc[idx,"select_on"]
        if not issubclass(type(v),str):
            err = "\n'select_on' should be string\n\n"
            raise ValueError(err)
        v = v.strip()

        # This will throw an error if not recognized
        ens.get_observable_function(obs_fcn=v)
        select_on_out.append(v)

    if len(set(select_on_out)) != 1:
        err = "\nselect_on must be the same for all conditions\n\n"
        raise ValueError(err)

    df["select_on"] = select_on_out

    # select_on_folded
    select_on_folded_out = []
    for idx in df.index:
        v = check_bool(value=df.loc[idx,"select_on_folded"],
                        variable_name="select_on_folded")
        select_on_folded_out.append(v)
    df["select_on_folded"] = select_on_folded_out        

    # temperature
    temperature_out = []
    for idx in df.index:
        v = check_float(value=df.loc[idx,"temperature"],
                        variable_name=f"temperature[{idx}]",
                        minimum_allowed=0,
                        minimum_inclusive=False)
        temperature_out.append(v)
    df["temperature"] = temperature_out

    # Construct a ligand dictionary
    ligand_dict = {}
    for c in ligands:

        # Make sure these are all floats
        for idx in df.index:
            v = check_float(value=df.loc[idx,c],
                            variable_name=f"ligand '{c}' concentration")
            df.loc[idx,c] = v
        ligand_dict[c] = df[c]

    return df, ligand_dict