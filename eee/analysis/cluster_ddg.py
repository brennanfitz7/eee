"""
Analyze a ddg file. 
"""

import eee
from eee._private.check.standard import check_int
from eee._private.check.standard import check_float

from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

from sklearn.cluster import KMeans
import numpy as np

# Choose whether to show (interactive) or close (non-interactive)
if eee._in_notebook:
    _plot_close_out = plt.show
else:
    _plot_close_out = plt.close

def _clean_up_axes(min_value,
                   max_value,
                   num_ticks=5,
                   scalar=0.05):
    """
    Clean up a set of axes, expanding min and max by scalar and then returning
    the new min/max with ticks. 
    """
    
    # min and max must be different
    if min_value == max_value:
        err = f"\nmin and max values are the same ({min_value})\n\n"
        raise ValueError(err)
    
    # force min_value to be less than max_value
    if min_value > max_value:
        _tmp = max_value
        max_value = min_value
        min_value = _tmp

    # Span 
    span = max_value - min_value

    # figure out how to expand min value
    min_value = min_value - span*scalar
    
    # figure out how to expand max value
    max_value = max_value + span*scalar
    
    # Get step. 
    step_size = span/(num_ticks-1)
    
    # Find ticks
    ticks = np.arange(num_ticks)*step_size + min_value
    
    return min_value, max_value, ticks
    

def cluster_ddg(ddg_df,
                max_num_clusters=20,
                elbow_cutoff=0.1,
                write_prefix="ddg-mut-effects",
                exclude_columns=None):
    """
    Cluster mutations based on their correlated effects on the different 
    species in the ennsemble.
    
    Parameters
    ----------
    ddg_df : str or pandas.DataFrame 
        dataframe with the energetic effects of mutations on species in an
        ensemble. 
    max_num_clusters : int, default=20
        when doing k-means clustering, try up to max_num_clusters clusters 
        when doing an elbow plot analysis. Must be positive int. 
    elbow_cutoff : float, default=0.1
        define breakpoint as point where adding a cluster only changes the 
        inertia by this fractional of the maximum inertia. Must be between 0 
        and 1.
    write_prefix : str, default="ddg-mut-effects"
        create graphical outputs with this file prefix. If None, no files are 
        written out.
    exclude_columns : list-like, optional
        drop these columns if they are in the dataframe. 
    """

    # Check inputs.

    input_df = eee.io.read_ddg(ddg_df)

    max_num_clusters = check_int(value=max_num_clusters,
                                 variable_name="max_num_clusters",
                                 minimum_allowed=1)
    
    elbow_cutoff = check_float(value=elbow_cutoff,
                               variable_name="elbow_cutoff",
                               minimum_allowed=0.0,
                               maximum_allowed=1.0)
    
    if write_prefix is not None:
        write_prefix = f"{write_prefix}"

    # Figure out which columns to drop. 
    drop_col = ["site","mut"]
    if exclude_columns is not None:
        exclude_columns = list(exclude_columns)
        drop_col.extend(exclude_columns)
    drop_col = list(set(drop_col).intersection(set(input_df.columns)))

    # Drop columns
    df = input_df.drop(columns=drop_col)
    
    # List of species in dataframe. 
    species = list(df.columns)
    species.sort()

    # ddG effects as an array. 
    data = np.array(df.loc[:,species])

    # Try from 1 to max_num_clusters clusters. 
    print("Finding appropriate number of clusters.",flush=True)
    num_clusters = np.arange(1,max_num_clusters + 1,dtype=int)

    # calculate k-means inertias. 
    inertias = []
    for i in num_clusters:
        kmeans = KMeans(n_clusters=i)
        kmeans.fit(data)
        inertias.append(kmeans.inertia_)

    # If only one cluster allowed, only one cluster
    if max_num_clusters == 1:
        found_number = 1

    # Find cluster using an elbow plot
    else:

        # Find fractional change in inertia with increased cluster number
        inertias = np.array(inertias)
        fx_changes = (inertias[1:-1] - inertias[2:])/inertias[1]

        # Appropriate cluster number from elbow plot
        mask = fx_changes > elbow_cutoff
        if sum(mask) == len(mask):
            found_number = max_num_clusters
        else:
            found_number = num_clusters[2:][mask][-1] - 1

    fig, ax = plt.subplots(1,figsize=(6,6))
    
    ax.plot(num_clusters, inertias, marker='o')
    ax.plot((found_number,found_number),(0,inertias[0]),"--",color='red')

    ax.set_title(r'$\Delta \Delta G$ elbow plot')
    ax.set_xlabel('number of clusters')
    ax.set_ylabel('inertia')

    if write_prefix is not None:
        fig.savefig(f"{write_prefix}_kmeans-elbow-plot.pdf")

    _plot_close_out()

    print("Clustering ddG values.",flush=True)
    
    kmeans = KMeans(n_clusters=found_number)
    kmeans.fit(data)
    
    fig, ax = plt.subplots(data.shape[1]-1,data.shape[1]-1,figsize=(10,10))

    for i in range(ax.shape[0]):
        for j in range(ax.shape[1]):
            ax[i,j].axis('off')

    for i in range(data.shape[1]):
        for j in range(i+1,data.shape[1]):

            ax[i,j-1].axis("on")
            ax[i,j-1].spines[['right', 'top']].set_visible(False)

            x = data[:,i]
            y = data[:,j]

            min_value = np.min([np.min(x),np.min(y)])
            max_value = np.min([np.max(x),np.max(y)])

            min_value, max_value, tick_positions = _clean_up_axes(min_value=min_value,
                                                                  max_value=max_value)

            draw_collection = ax[i,j-1].scatter(x, y, c=kmeans.labels_)

            ax[i,j-1].plot((min_value,max_value),(min_value,max_value),'--',color='gray')
            ax[i,j-1].set_xlabel(species[i] + r" $\Delta \Delta G$")
            ax[i,j-1].set_ylabel(species[j] + r" $\Delta \Delta G$")
            ax[i,j-1].set_xlim(min_value,max_value)
            ax[i,j-1].set_ylim(min_value,max_value)
            ax[i,j-1].set_xticks(tick_positions)
            ax[i,j-1].set_yticks(tick_positions)
            ax[i,j-1].set_aspect('equal', adjustable='box')

    ax[ax.shape[0]-1,0].set_axis_off()

    # Create legend color dict
    all_labels = list(np.unique(kmeans.labels_))
    all_labels.sort()
    color_dict = {}
    for k in all_labels:
        color_dict[k] = draw_collection.to_rgba(k)

    legend_entries = []
    for i in range(len(all_labels)):
        legend_entries.append(Line2D([0],[0],
                                     marker='o',
                                     color='w',
                                     label=f"cluster {all_labels[i]}",
                                     markerfacecolor=color_dict[i],
                                     markersize=10))

    ax[ax.shape[0]-1,0].legend(handles=legend_entries,
                               loc='center')

    plt.tight_layout()

    if write_prefix is not None:
        fig.savefig(f"{write_prefix}_clusters.pdf")

    _plot_close_out()

    columns = df.columns[:]

    df["cluster"] = kmeans.labels_

    for c in drop_col:
        df[c] = input_df[c]
        columns = columns.insert(0,c)

    columns = columns.insert(0,"cluster")

    df = df.loc[:,columns]
    
    return df
    
    