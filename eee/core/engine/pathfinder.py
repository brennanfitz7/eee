"""
Find accessible paths through sequence space starting at the wildtype genotype.
"""
from eee.core.genotype import Genotype

from eee._private.check.standard import check_bool
from eee._private.check.standard import check_int

from eee.analysis import get_num_genotypes

import numpy as np
from tqdm.auto import tqdm

def _traverse(current_genotype_idx,
              max_depth,
              current_depth,
              visited_dict,
              condition_fcn,
              not_exhaustive_path,
              gc,
              pbar):
    """
    Recursive function for exploring paths through sequence space. 
    current_genotype_idx is the current genotype, max_depth is how far we can 
    go, current_depth is the current distance from the start, visited_dict
    is a dictionary of visited genotypes (used to prevent visiting the same
    genotype by a different path), condition_fcn is a function that takes two
    values and returns True or False based on the comparison of those values 
    (_traverse uses the function as `condition_fcn(new_fitness,old_fitness)`),
    gc is a Genotype object, and pbar is a status bar. 
    """
    
    # We've reached the max depth; do not do anything
    if current_depth >= max_depth:
        return visited_dict
    
    # Consider all sites
    for site in gc.ddg_dict:
        
        # Can't mutate same site we have already mutated
        if site in gc.genotypes[current_genotype_idx].sites:
            continue
        
        # All possible mutations at this site
        for mutation in gc.ddg_dict[site]:
            
            # Construct mutant tuple and see if we've visited it before. 
            new_genotype = gc.genotypes[current_genotype_idx].mutations[:]
            new_genotype.append(mutation)
            new_genotype.sort()
            new_genotype = tuple(new_genotype)
        
            # Already tried this genotype
            if new_genotype in visited_dict and not_exhaustive_path:
                continue
            
            # Try to introduce the mutation
            new_idx = gc.conditional_mutate(index=current_genotype_idx,
                                            site=site,
                                            mutation=mutation,
                                            condition_fcn=condition_fcn)
            # Mutation accepted
            pbar.update(1)
            if new_idx >= 0:
                
                 # Record that we visited
                visited_dict[new_genotype] = True
                visited_dict = _traverse(current_genotype_idx=new_idx,
                                         max_depth=max_depth,
                                         current_depth=current_depth+1,
                                         visited_dict=visited_dict,
                                         condition_fcn=condition_fcn,
                                         not_exhaustive_path=not_exhaustive_path,
                                         gc=gc,
                                         pbar=pbar)
            else:
                visited_dict[new_genotype] = False

    return visited_dict


def pathfinder(gc,
               max_depth=1,
               allow_neutral=True,
               find_all_paths=True,
               output_file="genotypes.csv",
               return_output=False):
    """
    Recursively explore paths through sequence space implied by mutations stored
    in gc. 
    
    Parameters
    ----------
    gc : eee.core.Genotype
        initialized Genotype object with ensemble and ddg_df loaded
    max_depth : int, default=1s
            explore paths up to this length. 1 corresponds to all accessible
            single mutants, 2 to all double mutants, 3 to all triple, etc. 
    allow_neutral : bool, default=True
        allow mutations that have no effect on fitness
    find_all_paths : bool, default=True
        visit the same sequence from different starting points. This will find 
        all allowed paths to a given genotype. If False, ignore new paths that
        visit the same genotype. 
    output_file : str, default="genotypes.csv"
        return visited genotypes (with trajectory information) to this file
    return_output : bool, default=False
        return the dataframe describing the results. Otherwise, return None
    """

    if not issubclass(type(gc),Genotype):
        err = f"\ngc '{gc}' should be an instance of simulation.Genotype\n\n"
        raise ValueError(err)

    max_depth = check_int(value=max_depth,
                                variable_name="max_depth",
                                minimum_allowed=1,
                                minimum_inclusive=True)

    allow_neutral = check_bool(value=allow_neutral,
                               variable_name="allow_neutral")
    
    find_all_paths = check_bool(value=find_all_paths,
                                variable_name="find_all_paths")
    
    output_file = f"{output_file}"

    return_output = check_bool(value=return_output,
                               variable_name="return_output")
    
    if allow_neutral:
        condition_fcn = np.greater_equal
    else:
        condition_fcn = np.greater

    num_genotypes_per_shell = get_num_genotypes(gc.ddg_dict,max_depth=max_depth)
    total_calcs = np.sum(num_genotypes_per_shell) - 1

    pbar = tqdm(total=total_calcs)

    with pbar:

        print(f"Visiting genotypes up to {max_depth} mutations from wildtype.")
        if allow_neutral:
            print(f"Allowing adaptive and neutral steps")
        else:
            print(f"Allowing only adaptive steps")
        print("\nNOTE: status bar gives the maximum time if all genotypes are accessible.",flush=True)

        _ = _traverse(current_genotype_idx=0,
                    max_depth=max_depth,
                    current_depth=0,
                    visited_dict={},
                    condition_fcn=condition_fcn,
                    not_exhaustive_path=(not find_all_paths),
                    gc=gc,
                    pbar=pbar)
        
        # Indicate we completed the calculation
        pbar.n = total_calcs
        pbar.refresh()
        
    df = gc.df

    if not output_file is None:
        print("Writing output",flush=True)
        df.to_csv(output_file)

    out = None
    if return_output:
        out = df

    return out