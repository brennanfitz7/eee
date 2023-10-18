"""
Plot an energy landscape for an ensemble.
"""

from eee._private.check.ensemble import check_ensemble
from eee._private.check.plot import check_plot

import numpy as np
from matplotlib.lines import Line2D

def energy_diagram(values,
                   fig=None,
                   ax=None,
                   context_manager=None):
    """
    Plot an ensemble energy landscape.

    Parameters
    ----------
    ens : eee.core.Ensemble
        initialized ensemble with both observable and non-observable species
    fig : matplotlib.Figure, optional
        matplotlib Figure object. If this is specified, ax must be specified 
        as well. If fig is not specified and ax is, ignore fig argument. If
        neither fig nor ax are specified, generate them in this function. 
    ax : matplotib.Axes, optional
        matplotlib axis on which to generate the plot. If not specified, create
        a new axis object. 
    context_manager : object, optional
        do plot with this context manager. (Used for generating plots in 
        specific windows in a jupyte notebook). context_manager should be 
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
   
    #ens = check_ensemble(ens,check_obs=False)
    fig, ax, context_manager = check_plot(fig,ax,context_manager)

    with context_manager:

        all_dG0 = []
        mins_maxes = []
        for i, s in enumerate(values["species"]):

            dG0 = s["dG0"]
            if s["observable"]:
                color = "blue"
            else:
                color = "black"

            ax.plot((i,i+0.75),(dG0,dG0),'-',lw=2,color=color)

            if "ligands" in s:
                for lig in s["ligands"]:
                    if lig[0] == "" or lig[1] == 0:
                        continue

                    ax.arrow(x=i+0.375,
                             y=dG0,
                             dx=0,
                             dy=-lig[1],
                             length_includes_head=True,
                             head_width=0.02,
                             facecolor="black")
                    mins_maxes.append(dG0-lig[1])

            if not s["folded"]:
                ax.scatter(np.linspace(i,i+0.75,10),
                           np.ones(10)*dG0,
                           s=60,
                           marker="x",
                           color=color)

            all_dG0.append(dG0)
            mins_maxes.append(dG0)

        min_dG0 = np.min(mins_maxes)
        max_dG0 = np.max(mins_maxes)

        if max_dG0 - min_dG0 < 1.0:
            min_dG0 = min_dG0 - 0.5
            max_dG0 = max_dG0 + 0.5

        if min_dG0 > 0:
            min_dG0 = 0.95*min_dG0
        elif min_dG0 == 0:
            min_dG0 -= (max_dG0 - min_dG0)*0.05
        else:
            min_dG0 = 1.05*min_dG0

        if max_dG0 > 0:
            max_dG0 = 1.05*max_dG0
        elif max_dG0 == 0:
            max_dG0 += (max_dG0 - min_dG0)*0.05
        else:
            max_dG0 = 0.95*max_dG0

        span = max_dG0 - min_dG0

        ax.set_ylim(min_dG0,max_dG0)

        for i, s in enumerate(values["species"]):
            ax.text(i + 0.375,
                    all_dG0[i] + span*0.025,
                    s["name"],
                    size=14,
                    horizontalalignment="center")   

        ax.set_xlabel("species")
        ax.set_ylabel("free energy")
        ax.spines[['right', 'top','bottom']].set_visible(False)
        ax.set_xticks([])

        legend_elements = [Line2D([0],
                                  [0],
                                  color='blue',
                                  lw=2,
                                  label='observable'),
                           Line2D([0],
                                  [0],
                                  marker='x',
                                  color='black',
                                  label='unfolded',
                                  markersize=10)]

        ax.legend(handles=legend_elements)
        
    return fig, ax