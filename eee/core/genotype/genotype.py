"""
Class for keeping track of genotypes during an evolutionary simulation.
"""

from .single_genotype import SingleGenotype
from eee.io.read_ddg import read_ddg
from eee._private.check.ensemble import check_ensemble

import numpy as np
import pandas as pd

import os

class Genotype:
    """
    Class storing a collection of genotypes. It also allows mutation of one
    genotype into another.

    Key public attributes:

        + genotypes: all genotypes seen (as SingleGenotype instances).

        + trajectories: the trajectory taken by each genotype. Each entry
          is a list of steps taken to get to that genotype, where the steps are
          integer indexes pointing to self.genotype. A trajectory of [0,10,48]
          means the genotype started as self.genotypes[0], became
          self.genotypes[10], then finally self.genotypes[48].
          
        + mut_energies: mutational energies for the genotype.

        + fitnesses holds the absolute fitness of each genotype.
    """

    def __init__(self,ens,fitness_function,ddg_df,choice_function=None):
        """
        Initialize class. 

        Parameters
        ----------
        ens : Ensemble
            instance of Ensemble class
        fitness_function : function
            function that takes mutation energy and return fitness
        ddg_df : pandas.DataFrame
            dataframe holding the effects of each mutation on the energy of 
            each species. 
        choice_function : function, optional
            function that randomly selects elements from a list-like object 
            (i.e. np.random.choice). This can be passed in so we can use a 
            seeded and reproducible random number generator. 
        """
        
        # Fitness and ddg information
        self._ens = check_ensemble(ens,check_obs=True)
        self._fitness_function = fitness_function
        self._ddg_df = read_ddg(ddg_df)
    
        if choice_function is None:
            choice_function = np.random.choice
        self._choice_function = choice_function

        self._create_ddg_dict()

        # Sites and mutations for generating mutations
        self._possible_sites = list(self._ddg_dict.keys())
        self._mutations_at_sites = dict([(s,list(self._ddg_dict[s].keys()))
                                         for s in self._possible_sites])
        self._last_index = 0
        
        # Main public attributes of the class. 
        self._genotypes = {0:SingleGenotype(self._ens,self._ddg_dict)}
        self._trajectories = {0:[0]}
        self._mut_energies = {0:self._genotypes[0].mut_energy}

        wt_fitness = np.product(self._fitness_function(self._genotypes[0].mut_energy))
        self._fitnesses = {0:wt_fitness}

    def _create_ddg_dict(self):
        """
        Convert a ddg_df dataframe in to a dictionary of the form:
        
        ddg_dict["site"]["mut"] = np_array_with_mutant_effects_on_species
        """

        species = self._ens.species
        for s in species:
            if s not in self._ddg_df.columns:
                err = f"\nspecies {s} is not in ddg_df.\n\n"
                raise ValueError(err)
            
        self._ddg_dict = {}
        for i in self._ddg_df.index:
            row = self._ddg_df.loc[i,:]
            
            site = row["site"]
            mut = row["mut"]
            if site not in self._ddg_dict:
                self._ddg_dict[site] = {}
                
            mut_dict = {}
            for s in species:
                mut_dict[s] = row[s]
            
            self._ddg_dict[site][mut] = self._ens.mut_dict_to_array(mut_dict)


    def _add_genotype(self,new_genotype,prev_index):
        """
        Add a newly created genotype to the object. Updates _genotypes, 
        _mut_energies, _trajectories, and _fitnesses. Returns the key pointing
        to the new genotype.
        """

        # Get index for new genotype
        new_index = self._last_index + 1
        self._last_index += 1

        # Record the new genotype
        self._genotypes[new_index] = new_genotype

        # Record the trajectory taken to reach the current genotype
        new_trajectory = self._trajectories[prev_index][:]
        new_trajectory.append(new_index)
        self._trajectories[new_index] = new_trajectory
        
        self._mut_energies[new_index] = new_genotype.mut_energy

        # Record the fitness of this genotype
        new_fitness = np.product(self._fitness_function(new_genotype.mut_energy))
        self._fitnesses[new_index] = new_fitness

        return new_index

    def mutate(self,index,site=None,mutation=None):
        """
        Mutate the genotype with "index" to a new genotype, returning the 
        index of the new genotype.

        Parameters
        ----------
        index : int
            index of genotype to mutation. Must be a key in the
            self.genotypes dictionary.
        site : int, optional
            site to mutate. if not specified, choose one randomly. 
        mutation : str, optional
            mutation at the site to mutate. if not specified, choose one 
            randomly. 

        Returns
        -------
        new_key : int
            key of newly generated genotype
        """

        # Sanity check
        if index not in self._genotypes:
            err = f"\nindex ({index}) is not in genotypes\n\n"
            raise IndexError(err)
        
        # Create a new genotype and mutate
        new_genotype = self._genotypes[index].copy()

        # If not specified, randomly choose a site
        if site is None:
            site = self._choice_function(self._possible_sites)

        # If not specified, randomly choose a mutation
        if mutation is None:
            mutation = self._choice_function(self._mutations_at_sites[site])

        # Introduce mutation
        new_genotype.mutate(site,mutation)
        
        return self._add_genotype(new_genotype,index)


    def conditional_mutate(self,
                           index,
                           site,
                           mutation,
                           condition_fcn):
        """
        Introduce a mutation, making sure its fitness meets a defined condition.

        Parameters
        ----------
        index : int
            index of genotype to mutation. Must be a key in the
            self.genotypes dictionary.
        site : int, optional
            site to mutate. if not specified, choose one randomly. 
        mutation : str, optional
            mutation at the site to mutate. if not specified, choose one 
            randomly. 
        condition_fcn : function
            function that takes two floats as inputs an returns a bool. 
            np.equal would require the new genotype have the same fitness as 
            the parent genotype; np.greater would require the new genotype have
            improved fitness relative to the parent. 

        Returns
        -------
        new_index : int
            key for new genotype in self._genotypes. If the condition is not 
            met, do not store the new genotype and return -1. 
        """
        
        # Sanity check
        if index not in self._genotypes:
            err = f"\nindex ({index}) is not in genotypes\n\n"
            raise IndexError(err)
                           
        # Create a new genotype and mutate
        new_genotype = self._genotypes[index].copy()

        # Introduce mutation
        new_genotype.mutate(site,mutation)

        new_fitness = np.product(self._fitness_function(new_genotype.mut_energy))

        if condition_fcn(new_fitness,self._fitnesses[index]):
            return self._add_genotype(new_genotype,index)
        
        return -1
        


    def dump_to_csv(self,
                    filename,
                    keep_genotypes=None):
        """
        Dump genotypes into a csv file. This appends to the filename if the
        csv already exists. The dump removes the genotypes from the object. This
        is to allow long evolutionary simulations that could potentially take 
        a large amount of memory. If run without keep_genotypes specified, the
        resulting Genotype instance will be empty. 
        
        Parameters
        ----------
        filename : str
            csv file to write to
        keep_genotypes : list, optional
            list of genotypes to keep. These should be keys in self.genotypes    
        """

        # Get a dataframe describing the genotypes
        df = self.df

        # Create a genotype that only has genotypes to remove
        if keep_genotypes is not None:
            to_write = np.logical_not(df["genotype"].isin(keep_genotypes))
            df = df.loc[to_write,:]

         # Write out, either appending or creating a new file
        if os.path.isfile(filename):
            df.to_csv(filename,
                      mode="a",
                      header=False,
                      index=False)
        else:
            df.to_csv(filename,
                      mode="w",
                      header=True,
                      index=False)

        # Update genotypes, trajectories, mut_energies, and fitnesses. Default
        # is to drop all genotypes. Save those in keep_genotypes. 
        genotypes = {}
        trajectories = {}
        mut_energies = {}
        fitnesses = {}
        if keep_genotypes is not None:

            for g in keep_genotypes:
                genotypes[g] = self._genotypes[g]
                trajectories[g] = self._trajectories[g]
                mut_energies[g] = self._mut_energies[g]
                fitnesses[g] = self._fitnesses[g]
            
        self._genotypes = genotypes
        self._trajectories = trajectories
        self._mut_energies = mut_energies
        self._fitnesses = fitnesses

    def to_dict(self):
        """
        Return a json-able dictionary describing the ddg genotype parameters.
        """

        return {"ddg_df":self._ddg_df}

    @property
    def df(self):
        """
        Genotypes, trajectories, energies, and fitnesses as a pandas Dataframe.
        """

        genotypes = list(self._genotypes.keys())
        mutations = ["/".join(self._genotypes[g].mutations) for g in genotypes]
        num_mutations = [len(self._genotypes[g].mutations) for g in genotypes]
        accum_mutations = ["/".join(self._genotypes[g].mutations_accumulated)
                           for g in genotypes]
        num_accum_mut = [len(self._genotypes[g].mutations_accumulated)
                         for g in genotypes]

        parent = [self._trajectories[t][-2]
                  for t in list(self._trajectories.keys())[1:]]
        parent.insert(0,pd.NA)

        out = {"genotype":genotypes,
               "mutations":mutations,
               "num_mutations":num_mutations,
               "accum_mut":accum_mutations,
               "num_accum_mut":num_accum_mut,
               "parent":parent,
               "trajectory":[self._trajectories[t] for t in self._trajectories]}

        # Get mutation energies
        mut_energy_out = {}
        for name in self._ens.species:
            mut_energy_out[name] = []

        for i in self._genotypes:
            for j, name in enumerate(self._ens.species):
                mut_energy_out[name].append(self._mut_energies[i][j])

        for k in mut_energy_out:
            out[f"{k}_ddg"] = mut_energy_out[k]

        out["fitness"] = [self._fitnesses[f] for f in self._fitnesses]
        
        return pd.DataFrame(out)

    @property
    def wt_sequence(self):
        """
        Wildtype protein sequence (extracted from possible mutations).
        """

        wt = []
        for s in self._possible_sites:
            wt.append(self._mutations_at_sites[s][0][0])

        return "".join(wt)

    @property
    def genotypes(self):
        """
        Dictionary of genotypes.
        """

        return self._genotypes
    
    @property
    def trajectories(self):
        """
        Dictionary of trajectories that were taken on the way to each genotype.
        """
        return self._trajectories
    
    @property
    def mut_energies(self):
        """
        Dictionary of energies for each genotype.
        """
        return self._mut_energies

    @property
    def fitnesses(self):
        """
        Dictionary of absolute fitnesses.
        """
        return self._fitnesses
    
    @property
    def ddg_dict(self):
        """
        Dictionary with the energetic effects of mutations organized by site. 
        Has form :code:`ddg_dict[site][mut] = np_array_with_mutant_effects_on_species`,
        where site is an integer of sites from ddg_df, mut is a matched mutation
        at that site, and the array has the energetic effects of the mutation on 
        each species as ordered in ens.species. 
        """
        return self._ddg_dict



        


    