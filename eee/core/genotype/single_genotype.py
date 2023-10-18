"""
Class that holds an individual genotype.
"""

import numpy as np

import copy

class SingleGenotype:
    """
    Class to hold a genotype including the mutated sites, mutations, and the
    combined energetic effects of all mutations. Generally this will be 
    initialized by Genotype. 
    """

    def __init__(self,
                 ens,
                 ddg_dict,
                 sites=None,
                 mutations=None,
                 mutations_accumulated=None,
                 mut_energy=None):
        """
        Initialize instance.

        Parameters
        ----------
        ens : eee.core.Ensemble
            initialized ensemble
        ddg_dict : dict
            nested dictionary of mutational effects [site][mutation]
        sites : list, optional
            list of sites that already have mutations
        mutations : list, optional
            list of mutations to the genotype in an absolute Hamming distance 
            sense relative to wildtype. There is guaranteed to be a single 
            mutation at each site in this list. 
        mutations_accumulated : list, optional
            list of all mutations that have occurred to this genotype, in order, 
            over its history
        mut_energy : numpy.ndarray
            numpy array holding energetic effect of mutation(s) on each species.
            array should be in ensemble.species order. 
        """

        # These point to *instances* of these objects
        self._ens = ens
        self._ddg_dict = ddg_dict

        # NOTE: mutations and mutations_accumulated are shallow copies of lists
        # of strings. This is okay because all the lists do is declare which 
        # mutations are present in a genotype -- it's fine that all genotypes 
        # share individual string instances. It also makes the code way faster
        # because we do not need to copy the strings over and over. The sites
        # and mut_energy arrays are lists of numbers and numpy arrays and are
        # thus true copies. 

        # Sites is a new copy
        if sites is None:
            self._sites = []
        else:
            self._sites = copy.copy(sites)
        
        # Mutations is a new copy. 
        if mutations is None:
            self._mutations = []
        else:
            self._mutations = copy.copy(mutations)

        # Mutations is a new copy. 
        if mutations_accumulated is None:
            self._mutations_accumulated = []
        else:
            self._mutations_accumulated = copy.copy(mutations_accumulated)

        # mut_energy is a new copy. 
        if mut_energy is None:
            self._mut_energy = np.zeros(len(self._ens.species),dtype=float)
        else:
            self._mut_energy = mut_energy.copy()
            

    def copy(self):
        """
        Return a copy of the SingleGenotype instance. This will copy the ens and
        ddg_dict objects as references, but make new instances of the sites, 
        mutations, and mut_energy attributes. 
        """
        
        return SingleGenotype(self._ens,
                              self._ddg_dict,
                              sites=self._sites,
                              mutations=self._mutations,
                              mutations_accumulated=self._mutations_accumulated,
                              mut_energy=self._mut_energy)


    def mutate(self,site,mutation):
        """
        Introduce a mutation at a specific site. If this mutation has already
        occurred at this site, treat as a reversion. 
        """

        # Mutation to revert before adding new one
        prev_mut = None

        # If the site was already mutated, we need to mutate site back to wt
        # before mutating to new genotype
        if site in self._sites:

            # Get site of mutation in genotype
            idx = self._sites.index(site)

            # Get mutation to revert before adding new one
            prev_mut = self._mutations[idx]

            # Subtract the energetic effect of the previous mutation
            self._mut_energy -= self._ddg_dict[site][prev_mut]

            # Remove the old mutation
            self._sites.pop(idx)
            self._mutations.pop(idx)

        # If the mutation is the same as the previous mutation, treat as a
        # reversion back to the wildtype amino acid at this site. This means no
        # energetic effect of the mutation added to mut_energy and that the 
        # mutation and sites should be removed.
        if mutation != prev_mut:

            # Update energy of each species
            self._mut_energy += self._ddg_dict[site][mutation]

            # Update genotype with new site mutation and mutation made
            self._sites.append(site)
            self._mutations.append(mutation)

            # This was a mutation from the previous mutation state to the new 
            # state. Change mutation string to indicate
            final_mut_name = mutation[:]
            if prev_mut is not None:
                final_mut_name = f"{prev_mut[-1]}{mutation[1:]}"
        
        else:
            # This was a reversion, update the mutation string to indicate this
            final_mut_name = f"{mutation[-1]}{mutation[1:-1]}{mutation[0]}"
        
        # Record mutation that occurred (new site, change at existing site, 
        # reversion at site)
        self._mutations_accumulated.append(final_mut_name)


    @property
    def sites(self):
        """
        List of sites that are mutated relative to wildtype. Index synced with
        mutations.
        """
        return self._sites

    @property
    def mutations(self):
        """
        List of mutations in this genotype relative to wildtype. Index synced
        with sites.
        """
        return self._mutations
    
    @property
    def mutations_accumulated(self):
        """
        List of all mutations in the order they occurred.
        """
        return self._mutations_accumulated

    @property
    def mut_energy(self):
        """
        Current energetic effect of all mutations on the species in the 
        ensemble. This is a numpy array, with species order corresponding to 
        the ensemble.species order. 
        """
        return self._mut_energy
    
