"""
Determine the number of genotypes in each shell of mutational steps away 
from a central genotype. 
"""

from eee._private.check.standard import check_int

import numpy as np

import itertools
import math

def get_num_genotypes(ddg_dict,max_depth):
    """
    Determine the number of genotypes in each shell of mutational steps away 
    from a central genotype. 

    Parameters
    ----------
    ddg_dict : dict
        nested dictionary of mutational effects [site][mutation]
    max_depth : int
        number of steps away from wildtype to proceed

    Returns
    -------
    num_genotypes_in_shell : numpy.ndarray
        number of genotypes from 0 to max_depth genotypes away from wildtype. 
        The first shell has 1 genotype (wildtype); subsequent shells have how 
        many combinations of 1, 2, 3, etc. mutations are possible. 
    """

    if not issubclass(type(ddg_dict),dict):
        err = f"\nddg_dict '{ddg_dict}' should be a dictionary\n\n"
        raise ValueError(err)
    
    max_depth = check_int(max_depth,
                          variable_name=max_depth,
                          minimum_allowed=0)
    
    # Get sites to mutate
    all_sites = list(ddg_dict.keys())
    num_sites = len(all_sites)

    # Can't have more than num_sites mutations at once...
    if max_depth > num_sites:
        max_depth = num_sites


    # Figure out number of mutations at each site
    site_lengths = []
    for s in all_sites:
        site_lengths.append(len(ddg_dict[s]))
    site_lengths = np.array(site_lengths)

    num_genotypes_in_shell = [1]

    # Same number of mutations at each site, single loop 
    if len(np.unique(site_lengths)) == 1:
        N = site_lengths[0]

        for L in range(1,max_depth+1):
            num_genotypes = math.comb(num_sites,L)*(N**L)
            num_genotypes_in_shell.append(num_genotypes)


    # Different sites have different numbers of mutations -- brute force it
    else:

        for L in range(1,max_depth+1):
            num_genotypes = 0
            for sites in itertools.combinations(all_sites,L):
                num_genotypes += np.prod(site_lengths[np.array(sites)-1])
            num_genotypes_in_shell.append(num_genotypes)

    return np.array(num_genotypes_in_shell)