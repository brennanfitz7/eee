"""
Functions for summarizing the result of a Wright-Fisher simulation. 
"""

# Need to think through how to have unified plot rc (or even if we should) 
# across package. this should bring in matplotlib rc settings from basic plot 
# interface. 
from eee.analysis.plots import plt

from eee.calcs.read_json import read_json
from eee._private.interface import MockContextManager

import pandas as pd
from tqdm.auto import tqdm
import numpy as np

import os
import glob
import pickle


def plot_fitness(fit_df,
                 color='red',
                 label=None,
                 fig=None,
                 ax=None):
    """
    Plot population fitness versus generation.

    Parameters
    ----------
    fit_df : pandas.DataFrame   
        dataframe with "mu" and "std" columns corresponding to mean and standard
        deviation of fitness over generations
    color : str, default="red"
        color for mean line
    label : str, optional
        label (for legend assignment)
    fig : matplotlib.Figure, optional
        if specified, plot using this fig
    ax : matplotlib.Axis, optional
        if specified, plot using this axis

    Returns
    -------
    fig, ax : matplotlib.Figure, matplotlib.Axis
        figure and axis from plot. if these were passed in, they are passed back
        out. if not passed in, these are generated locally. 
    """


    if fig is None:
        fig, ax = plt.subplots(1,figsize=(8,5))

    x = np.arange(len(fit_df["mu"]))

    top_of_area = fit_df["mu"] + fit_df["std"]
    top_of_area[top_of_area > 1.0] = 1.0

    bottom_of_area = fit_df["mu"] - fit_df["std"]
    bottom_of_area[bottom_of_area < 0.0] = 0.0

    ax.fill_between(x,top_of_area,bottom_of_area,color='0.6')
    ax.plot(fit_df["mu"],lw=3,zorder=20,color=color,label=label)

    ax.set_ylim(np.min(fit_df["mu"])*0.95,1.05)
    ax.set_xlabel("generation")
    ax.set_ylabel("fitness")
    
    return fig, ax

def plot_genotypes_vs_generations(genotypes,
                                  traj_df,
                                  plot_top=None,
                                  fig=None,
                                  ax=None):
    """
    Plot genotype frequency versus time. 

    Parameters
    ----------
    genotypes : pandas.DataFrame
        genotype output from simulation
    traj_df : pandas.DataFrame
        trajectory dataframe where columns are genotypes and rows are generations.
        values are freq of that genotype at that generation
    plot_top : int, optional
        plot only the first plot_top most common genotypes
    fig : matplotlib.Figure, optional
        if specified, plot using this fig
    ax : matplotlib.Axis, optional
        if specified, plot using this axis

    Returns
    -------
    fig, ax : matplotlib.Figure, matplotlib.Axis
        figure and axis from plot. if these were passed in, they are passed back
        out. if not passed in, these are generated locally. 
    """

    if fig is None:
        fig, ax = plt.subplots(1,figsize=(6,6))

    to_plot = list(traj_df.columns)
    to_pool = ["other"]
    if "other" in to_plot:
        to_plot.remove("other")

    if plot_top is not None:
        to_plot = to_plot[:plot_top]
        to_pool.extend(to_plot[plot_top:])

    all_others = np.sum(traj_df.loc[:,to_pool],axis=1)
    
    for t in to_plot:

        mutations = genotypes.loc[int(t),"mutations"]
        if pd.isnull(mutations):
            mutations = "wt"
        fitness = genotypes.loc[int(t),"fitness"]

        ax.plot(traj_df[t],lw=2,label=f"{mutations}: {np.log10(fitness):.5f}")

    ax.plot(all_others,lw=2,color="gray",label="all others")

    ax.legend()
    ax.set_xlim(0,len(traj_df))
    ax.set_ylim(-0.05,1.05)
    ax.set_xlabel("generation")
    ax.set_ylabel("fraction of population")
    
    return fig, ax

def load_simulation_results(sim_directory):
    """
    Load the contents of a simulation into a dictionary. 
    
    Parameters
    ----------
    sim_directory : str
        directory holding output from a simulation
    
    Returns
    -------
    results : dict
        dictionary with keys holding Simulation object, keywords passed to 
        start the run, the genotypes dataframe, and a list of (sorted) pickle
        files holding the populations at each generation. 
    """
    
    so, kwargs = read_json(os.path.join(sim_directory,
                                        "input",
                                        "simulation.json"),
                                        use_stored_seed=True)

    if so.calc_type != "wf_sim":
        err = f"\nCalculation is type '{so.calc_type}', not 'wf_sim'.\n\n"
        raise ValueError(err)
    
    # Read dataframe
    df_file = os.path.join(sim_directory,f"{kwargs['write_prefix']}_genotypes.csv")
    genotypes = pd.read_csv(df_file)
    genotypes = genotypes.set_index("genotype")

    # Get list of pickle files
    pickle_pattern = os.path.join(sim_directory,f"{kwargs['write_prefix']}*.pickle")
    pickle_files = glob.glob(pickle_pattern)
    pickle_files.sort()
    
    results = {"so":so,
               "kwargs":kwargs,
               "genotypes":genotypes,
               "pickle_files":pickle_files}
    
    return results

def fitness_and_variance(genotypes,
                         pickle_files,
                         verbose=True):
    """
    Given the output of a wf_sim run, calculate the mean and standard deviation
    of the population fitness over time.
    
    Parameters
    ----------
    genotypes : pandas.DataFrame
        genotypes dataframe created by load_directory
    pickle_files : list
        list of pickle files to load (in order)
    verbose : bool, default=True
        print outputs and status bars
    
    Returns
    -------
    df : pandas.DataFrame
        dataframe with mean and standard deviation of the fitness
    """

    if verbose:
        print("Calculating fitness over time")
        print("Reading pickle files",flush=True)
        pbar = tqdm(total=len(pickle_files))
    else:
        pbar = MockContextManager()
        

    mean_fitness = []
    sd_fitness = []
    num_generations = 0 
    with pbar:

        for pickle_file in pickle_files:
    
            # Load this block of generations
            with open(pickle_file,'rb') as f:
                generations = pickle.load(f)

            # Get population size from first entry
            if num_generations == 0:
                population_size = np.sum(list(generations[0].values()))

            # Go through generations
            for generation in generations:
                
                # Create a vector holding the fitness of every genotype in this
                # generation
                last_counter = 0
                fitness_vector = np.zeros(population_size,dtype=float)
                for genotype in generation:
                    
                    f = genotypes.loc[int(genotype),"fitness"]
                    new_counter = last_counter + generation[genotype]
                    fitness_vector[last_counter:new_counter] = f
                    last_counter = new_counter
                
                # Record mean and sd of fitness vector
                mean_fitness.append(np.mean(fitness_vector))
                sd_fitness.append(np.std(fitness_vector))
                
                num_generations += 1
            
            pbar.update()

    return pd.DataFrame({"mu":mean_fitness,
                         "std":sd_fitness})
            

def extract_key_trajectories(pickle_files,
                             key_genotype_cutoff=1e-3,
                             verbose=True):
    """
    Extract frequency of genotypes seen in at least cutoff of the individuals 
    seen over the entire simulation. 

    Parameters
    ----------
    pickle_files : list
        list of pickle files (in order)
    key_genotype_cutoff : float, default=1e-3
        how to filter genotypes. if cutoff is 0.001, this will return the 
        genotypes that account for 0.999 of all individuals seen over the 
        entire simulation
    verbose : bool, default=True
        print information and status bars
    
    Returns
    -------
    out_df : pandas.DataFrame
        dataframe where columns are genotypes (can be looked up in the
        simulation genotypes dataframe) and rows are generations. entries are 
        frequency at that generation. the "other" columns accounts for genotypes
        discarded. 
    """

    if verbose:
        print("Extracting key trajectories")
        print("Counting the number of times each genotype was seen",flush=True)
        pbar = tqdm(total=len(pickle_files))
    else:
        pbar = MockContextManager()

    # Count how many times each genotype is seen over the entire trajectory
    num_generations = 0
    times_seen = {}
    with pbar:

        # Go through all generations
        for pickle_file in pickle_files:
            with open(pickle_file,'rb') as f:
                generations = pickle.load(f)

            # For each generation
            for generation in generations:

                # Get population size
                if num_generations == 0:
                    population_size = np.sum(list(generation.values()))

                # Build up times_seen count for how often we say a given 
                # genotype.
                for genotype in generation:
                    if genotype not in times_seen:
                        times_seen[genotype] = 0
                    times_seen[genotype] += generation[genotype]

                num_generations += 1
        
            pbar.update()

    # Record the total number of generations seen
    total_num_generations = num_generations

    if verbose:
        percentile = 100*(1 - key_genotype_cutoff)
        print(f"Identifying the set of genotypes that make up {percentile:.5f}% of individuals over time",flush=True)

    # Build cumulative distribution of genootypes
    all_counts = np.array(list(times_seen.values()))
    all_counts.sort()
    cumulative = np.cumsum(all_counts)/np.sum(all_counts)

    # Find count cutoff that corresponds to cutoff
    fx = np.sum(cumulative > key_genotype_cutoff)/len(cumulative)
    count_cutoff = np.min(all_counts[cumulative > key_genotype_cutoff])

    # Keep only genotypes seen at least count_cutoff times
    filtered_genotypes = {}
    for t in times_seen:
        if times_seen[t] >= count_cutoff:
            filtered_genotypes[t] = times_seen[t]

    if verbose:
        print(f" + Taking {len(filtered_genotypes)} ({100*fx:.3f}%) of genotypes")
        print(f" + Corresponds to genotypes seen >= {count_cutoff} times",flush=True)
    
        print("Setting up structures to hold genotype histories",flush=True)
        pbar = tqdm(total=len(filtered_genotypes)*2)
    else:
        pbar = MockContextManager()

    to_sort = []
    with pbar:
        for f in filtered_genotypes:
            to_sort.append((filtered_genotypes[f],f))
            pbar.update()
    
        to_sort.sort()
        to_sort = to_sort[::-1]

        genotype_traj = {}
        for t in to_sort:
            genotype_traj[t[1]] = np.zeros(num_generations,dtype=int)
            pbar.update()

    if verbose:
        print("Extracting genotype histories",flush=True)
        pbar = tqdm(total=len(pickle_files))
    else:
        pbar = MockContextManager()

    # Construct a dictionary to hold frequency of each genotype over generations
    to_df = dict([(f"{g}",np.zeros(total_num_generations,dtype=float))
                  for g in filtered_genotypes])
    with pbar:

        # Go through simulation pickles
        generation_counter = 0
        for pickle_file in pickle_files:

            with open(pickle_file,'rb') as f:
                generations = pickle.load(f)

            # For each generation
            for generation in generations:

                # For each genotype
                for genotype in generation:

                    # If the genotype is one we care about, record its frequency
                    if genotype in filtered_genotypes:
                        freq = generation[genotype]/population_size
                        to_df[f"{genotype}"][generation_counter] = freq

                generation_counter += 1
            
            pbar.update()

    # Create final dataframe
    out_df = pd.DataFrame(to_df)
    out_df["other"] = 1 - np.sum(out_df,axis=1)

    return out_df

def summarize_simulation(sim_directory,
                         key_genotype_cutoff=1e-3,
                         verbose=True,
                         generate_plots=True,
                         delete_pickle_files=False):
    """
    Summarize the results of a simulation, optionally deleting the raw pickle
    files.
    
    Parameters
    ----------
    sim_directory : str 
        directory holding simulation results
    key_genotype_cutoff : float
        how to filter genotypes. if cutoff is 0.001, this will return the 
        genotypes that account for 0.999 of all individuals seen over the 
        entire simulation
    verbose : bool, default=True
        print output and status bars
    generate_plots : bool, default=True
        generate plots and ave as pdf files in sim_directory
    delete_pickle_files : bool, default=False
        delete the pickle files in the directory
    
    Returns
    -------
    sim : dict
        simulation results with keys "so" (simulation object), "kwargs"
        (keyword arguments used on so.run), "genotypes" (dataframe with all 
        genotypes seen over simulation) and "pickle_files" (pickle files 
        written out by the simulation)
    fitness_results : pandas.DataFrame
        dataframe holding mean and variance of population fitness over the 
        simulation
    key_traj : pandas.DataFrame
        frequencies of genotypes that make up 1-cutoff of all genotypes seen
        as a function of time

    Notes
    -----
    This function also writes avg_fitness.csv and key_traj.csv to the simulation
    directory. If generate_plots == True, plots are written as pdf files to that
    directory as well. 
    """

    sim = load_simulation_results(sim_directory)
    fitness_results = None
    key_traj = None

    # Get overall fitness versus time
    avg_fitness_csv = os.path.join(sim_directory,"avg_fitness.csv")
    if os.path.exists(avg_fitness_csv):
        fitness_results = pd.read_csv(avg_fitness_csv)

    else:
        if len(sim["pickle_files"]) > 0:

            # We have pickle files. Load results and generate plots if requested. 
            fitness_results = fitness_and_variance(genotypes=sim["genotypes"],
                                                   pickle_files=sim["pickle_files"],
                                                   verbose=verbose)
            fitness_results.to_csv(avg_fitness_csv,index=False)
            
            if generate_plots:
                fig, _ = plot_fitness(fit_df=fitness_results)
                fig.savefig(os.path.join(sim_directory,"avg_fitness.pdf"))
                plt.close()
    
    # Get key trajectories
    key_traj_csv = os.path.join(sim_directory,"key_traj.csv")
    if os.path.exists(key_traj_csv):
        key_traj = pd.read_csv(key_traj_csv)

    else:
        if len(sim["pickle_files"]) > 0:
            
            # We have pickle files. Load results and generate plots if requested. 
            key_traj = extract_key_trajectories(pickle_files=sim["pickle_files"],
                                                key_genotype_cutoff=key_genotype_cutoff,
                                                verbose=verbose)
            key_traj.to_csv(os.path.join(sim_directory,"key_traj.csv"),
                            index=False)
            
            if generate_plots:
                fig, _ = plot_genotypes_vs_generations(genotypes=sim["genotypes"],
                                                       traj_df=key_traj)
                fig.savefig(os.path.join(sim_directory,"key_traj.pdf"))
                plt.close()
        

    # Delete pickle files, if requested
    if delete_pickle_files:
        
        if len(sim["pickle_files"]) > 0:

            # Read last pickle file
            with open(sim["pickle_files"][-1],'rb') as f:
                generations = pickle.load(f)
            
            # Write last generation from last pickle file to it's own pickle file
            # to allow extension of trajectories
            last_generation = [generations[-1]]
            with open(os.path.join(sim_directory,"last-generation.pickle"),"wb") as f:
                pickle.dump(last_generation,f)

        # Remove original pickle files
        for p in sim["pickle_files"]:
            os.remove(p)

    return sim, fitness_results, key_traj