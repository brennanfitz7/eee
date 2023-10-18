"""
Validate an ensemble is ready for a calculation. 
"""

from eee.core import Ensemble

def check_ensemble(ens,check_obs=False):
    """
    Validate an ensemble.
    
    Parameters
    ----------
    ens : eee.core.Ensemble
        ensemble to validate
    check_obs : bool, default=False
        make sure the observable function can be used. This will raise a 
        ValueError if the ensemble is not properly initialized.
    
    Returns
    -------
    ens : eee.core.Ensemble
        pointer to the ensmeble passed into the function. 
    """

    if not issubclass(type(ens),Ensemble):
        err = "ens should be an instance of the Ensemble class"
        raise ValueError(f"\n{err}\n\n")
    
    # This will throw an error if the ensemble is improperly initialized 
    # for an observable calculation (too few species, no observable species,
    # etc.)
    if check_obs:
        ens.get_obs()

    return ens