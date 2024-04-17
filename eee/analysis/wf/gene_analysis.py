"""
Scripts for analyzing the results of Wright-Fisher evolutionary simulations. 
"""

import eee

import numpy as np
import pandas as pd
from tqdm.auto import tqdm

from io import StringIO

def read_genotypes_file(genotypes_file,
                        only_genotypes=None):
    """
    Read a genotypes file.
    
    Parameters
    ----------
    genotypes_file : str
        path to genotypes file written out in a Wright-Fisher simulation
    only_genotypes : iterable, optional
        only get the genotypes in a list. should be a list of integers. Used to
        avoid loading a huge genotypes file into memory
    
    Returns
    -------
    genotypes : pandas.DataFrame
        genotypes as a dataframe
    """
    
    print(f"Reading genotypes file {genotypes_file}",flush=True)

    if only_genotypes is None:
        genotypes = eee.io.read_dataframe(genotypes_file)
        return genotypes
    
    only_genotypes = set(only_genotypes)

    # Count the number of lines in the file for a status bar
    def _make_gen(reader):
        b = reader(1024 * 1024)
        while b:
            yield b
            b = reader(1024*1024)

    f = open(genotypes_file, 'rb')
    f_gen = _make_gen(f.raw.read)
    num_lines = sum(buf.count(b'\n') for buf in f_gen) - 1

    pbar = tqdm(total=num_lines)
    with pbar:

        to_read = []
        with open(genotypes_file) as f:
            for line in f:
                if len(to_read) == 0:
                    to_read.append(line)
                    continue
                
                col = line.split(",")
                if int(col[0]) in only_genotypes:
                    to_read.append(line)
            
                pbar.update(n=1)

    genotypes = pd.read_csv(StringIO("".join(to_read)),sep=",")
    return genotypes

def get_genotype_frequencies(generations,cutoff=50):
    """
    Extract genotype frequencies over simulation. 
    
    Parameters
    ----------
    generations : list
        list of dicts holding results returned by a Wright-Fisher simulation
    cutoff : int, default=50
        drop genotypes seen fewer than cutoff times over the simulation
    
    Returns
    -------
    df : pandas.DataFrame
        dataframe with genotype frequencies over the course of the simulation. 
        columns are genotypes; rows are generations. 
    """

    max_counts = {}
    population_size = None
    for g in tqdm(generations):

        if population_size is None:
            population_size = 0
            for k in g:
                population_size += g[k]
                
        for k in g:
            if k not in list(max_counts.keys()):
                max_counts[k] = g[k]
                
            if g[k] > max_counts[k]:
                max_counts[k] = g[k]
                
    keys_to_keep = []
    for k in max_counts:
        if max_counts[k] > cutoff:
            keys_to_keep.append(k)
    
    out_dict = {}
    for k in keys_to_keep:
        out_dict[k] = np.zeros(len(generations),dtype=float)
        for i in range(len(generations)):
            try:
                out_dict[k][i] = generations[i][k]/population_size
            except KeyError:
                pass
        
    return pd.DataFrame(out_dict)
    

def get_most_common(some_dict):
    """
    Return key/value pair where value is highest. 
    """
    
    keys = list(some_dict.keys())
    values = [some_dict[k] for k in keys]
    to_sort = list(zip(values,keys))
    to_sort.sort()
    
    return to_sort[-1][1], to_sort[-1][0]
