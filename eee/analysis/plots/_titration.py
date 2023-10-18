"""
Plot a titration of ensemble observable versus ligand concentration. 
"""

from eee._private.check.ensemble import check_ensemble
from eee._private.check.plot import check_plot
from eee._private.check.eee import check_ligand_dict


def titration(ens,
              ligand_dict,
              fig=None,
              ax=None,
              context_manager=None):
    """
    Plot a titration of fx_obs versus concentration specified in ligand_dict.

    Parameters
    ----------
    ens : eee.core.Ensemble
        initialized ensemble with both observable and non-observable species
    ligand_dict : dict
        dictionary of chemical potentials. keys are the names of chemical
        potentials. Values are floats or arrays of floats. Any arrays 
        specified must have the same length. If a chemical potential is not
        specified in the dictionary, its value is set to 0. 
    fig : matplotlib.Figure, optional
        matplotlib Figure object. If this is specified, ax must be specified 
        as well. If fig is not specified and ax is, ignore fig argument. If
        neither fig nor ax are specified, generate them in this function. 
    ax : matplotib.Axes, optional
        matplotlib axis on which to generate the plot. If not specified, create
        a new axis object. 
    context_manager : object, optional
        do plot with this context manager. (Used for generating plots in 
        specific windows in a jupyter notebook). context_manager should be 
        initialized object with __enter__ and __exit__ methods. 

    Returns
    -------
    fig : matplotlib.Figure or None
        figure used for plot. If fig passed in as argument, this is returned. 
        If ax passed in alone, return None. If neither fig nor ax passed in, 
        return generated fig. 
    ax : matplotlib.Axes 
        axes used for plot. If ax was passed in as an argument, return this. 
        If ax was not passed in, return generated ax. 
    """

    ens = check_ensemble(ens,check_obs=True)
    ligand_dict, _ = check_ligand_dict(ligand_dict)
    fig, ax, context_manager = check_plot(fig=fig,
                                          ax=ax,
                                          context_manager=context_manager)

    with context_manager:
    
        ligand_key = list(ligand_dict.keys())[0]
        df = ens.get_obs(ligand_dict=ligand_dict)

        ax.plot(df[ligand_key],
                df["fx_obs"],
                color="blue",
                lw=2)

        ax.set_ylim(-0.05,1.05)
        ax.set_ylabel("fx observable")
        ax.set_xlabel(f"{ligand_key} chemical potential")
        ax.spines[['right', 'top']].set_visible(False)
    
    return fig, ax