import pytest

from eee.core.engine.pathfinder import _traverse
from eee.core.engine.pathfinder import pathfinder
from eee._private.interface import MockContextManager

import pandas as pd
import numpy as np

import copy 
import os
import glob

def test__traverse():

    # Line ensures test is seen by completeness_crawler
    assert True
    
    class _FakeSingleGenotype:
        """
        Fake class that exposes what we need to test _traverse.
        """

        def __init__(self,sites,mutations,fitness):
            self.sites = sites
            self.mutations= mutations
            self.fitness = fitness

    class _FakeGenotype:
        """
        Fake class that exposes what we need to test _traverse. mut_dict is a
        bit different -- this stores the fitness effect of each mutation at 
        each site. Mutational effects are entirely additive.
        """

        def __init__(self,ddg_dict):
            
            self.genotypes = {0:_FakeSingleGenotype(sites=[],
                                                    mutations=[],
                                                    fitness=0)}
            self.ddg_dict = ddg_dict

        def conditional_mutate(self,index,site,mutation,condition_fcn):

            F1 = self.genotypes[index].fitness
            F2 = self.genotypes[index].fitness + self.ddg_dict[site][mutation]
        
            if condition_fcn(F2,F1):
                new_idx = np.max(list(self.genotypes.keys())) + 1

                sites = self.genotypes[index].sites[:]
                sites.append(site)

                mutations = self.genotypes[index].mutations[:]
                mutations.append(mutation)

                self.genotypes[new_idx] = _FakeSingleGenotype(sites=sites,
                                                              mutations=mutations,
                                                              fitness=F2)
                return new_idx
            
            return -1

    def _always_allow(v1,v2):
        return True
    
    def _never_allow(v1,v2):
        return False

    def _build_genotype_tuple(genotypes):
        """
        List of tuples, one for each genotype, sorted.
        """
        
        out = []
        for g in genotypes:
            out.append(genotypes[g].mutations[:])
            out[-1].sort()
            out[-1] = tuple(out[-1])

        return out

    # Fake progress bar
    pbar = MockContextManager()

    # --------------------- allow any step ---------------------------

    ddg_dict = {1:{"A1V":1,"A1P":-10},
                2:{"A2C":0}}

    gc = _FakeGenotype(ddg_dict=ddg_dict)
        
    visited_dict = _traverse(current_genotype_idx=0,
                             max_depth=2,
                             current_depth=0,
                             visited_dict={},
                             condition_fcn=_always_allow,
                             not_exhaustive_path=True,
                             gc=gc,
                             pbar=pbar)
    
    assert len(visited_dict) == 5
    expected_visits = [("A1V",),("A1P",),("A2C",),("A1V","A2C",),("A1P","A2C")]
    for k in expected_visits:
        assert k in visited_dict

    assert len(gc.genotypes) == 6
    sorted_genotypes = _build_genotype_tuple(gc.genotypes)
    expected = [(),
                ('A1V',),('A1V','A2C'),
                ('A1P',),('A1P','A2C'),
                ('A2C',)]
    for i in range(len(sorted_genotypes)):
        assert np.array_equal(sorted_genotypes[i],expected[i])

    # --------------------- do not allow any step ---------------------------

    gc = _FakeGenotype(ddg_dict=ddg_dict)
    
    visited_dict = _traverse(current_genotype_idx=0,
                             max_depth=2,
                             current_depth=0,
                             visited_dict={},
                             condition_fcn=_never_allow,
                             not_exhaustive_path=True,
                             gc=gc,
                             pbar=pbar)
    
    assert len(visited_dict) == 3
    expected_visits = [("A1V",),("A1P",),("A2C",)]
    for k in expected_visits:
        assert k in visited_dict

    assert len(gc.genotypes) == 1
    sorted_genotypes = _build_genotype_tuple(gc.genotypes)
    expected = [()]
    for i in range(len(sorted_genotypes)):
        assert np.array_equal(sorted_genotypes[i],expected[i])


    # --------------------- use real fitness values ---------------------------

    gc = _FakeGenotype(ddg_dict=ddg_dict)
    
    visited_dict = _traverse(current_genotype_idx=0,
                             max_depth=2,
                             current_depth=0,
                             visited_dict={},
                             condition_fcn=np.greater_equal,
                             not_exhaustive_path=True,
                             gc=gc,
                             pbar=pbar)
    
    assert len(visited_dict) == 5
    expected_visits = [("A1V",),("A1P",),("A2C",),("A1V","A2C",),("A1P","A2C")]
    for k in expected_visits:
        assert k in visited_dict

    assert len(gc.genotypes) == 4
    sorted_genotypes = _build_genotype_tuple(gc.genotypes)
    expected = [(),
                ('A1V',),('A1V','A2C'),
                ('A2C',)]
    for i in range(len(sorted_genotypes)):
        assert np.array_equal(sorted_genotypes[i],expected[i])


    # --------------------- use real fitness values ---------------------------

    gc = _FakeGenotype(ddg_dict=ddg_dict)
    
    visited_dict = _traverse(current_genotype_idx=0,
                             max_depth=2,
                             current_depth=0,
                             visited_dict={},
                             condition_fcn=np.greater,
                             not_exhaustive_path=True,
                             gc=gc,
                             pbar=pbar)
    
    assert len(visited_dict) == 4
    expected_visits = [("A1V",),("A1P",),("A2C",),("A1V","A2C",)]
    for k in expected_visits:
        assert k in visited_dict

    assert len(gc.genotypes) == 2
    sorted_genotypes = _build_genotype_tuple(gc.genotypes)
    expected = [(),('A1V',)]
    for i in range(len(sorted_genotypes)):
        assert np.array_equal(sorted_genotypes[i],expected[i])

    # --------------------- Switch to exhaustive paths ---------------------------

    gc = _FakeGenotype(ddg_dict=ddg_dict)
    
    visited_dict = _traverse(current_genotype_idx=0,
                             max_depth=2,
                             current_depth=0,
                             visited_dict={},
                             condition_fcn=np.greater_equal,
                             not_exhaustive_path=False,
                             gc=gc,
                             pbar=pbar)
    
    assert len(visited_dict) == 5
    expected_visits = [("A1V",),("A1P",),("A2C",),("A1V","A2C",),("A1P","A2C")]
    for k in expected_visits:
        assert k in visited_dict

    assert len(gc.genotypes) == 5
    sorted_genotypes = _build_genotype_tuple(gc.genotypes)
    expected = [(),('A1V',),('A1V','A2C'),('A2C',),('A1V','A2C')]
    for i in range(len(sorted_genotypes)):
        assert np.array_equal(sorted_genotypes[i],expected[i])


def test_pathfinder(ens_with_fitness_two_site,variable_types,tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)
    
    gc = copy.deepcopy(ens_with_fitness_two_site["gc"])
    
    # Make sure it wrote out results
    out1 = pathfinder(gc=gc,
                      max_depth=1,
                      allow_neutral=True,
                      find_all_paths=False,
                      output_file="genotypes1.csv",
                      return_output=True)
    assert os.path.exists("genotypes1.csv")
    assert issubclass(type(out1),pd.DataFrame)
    
    # Max length increased
    gc = copy.deepcopy(ens_with_fitness_two_site["gc"])
    out2 = pathfinder(gc=gc,
                      max_depth=2,
                      allow_neutral=True,
                      find_all_paths=False,
                      output_file="genotypes2.csv",
                      return_output=True)
    assert os.path.exists("genotypes2.csv")
    assert issubclass(type(out2),pd.DataFrame)
    assert len(out1) < len(out2)
    
    # No neutral allowed
    gc = copy.deepcopy(ens_with_fitness_two_site["gc"])
    out3 = pathfinder(gc=gc,
                      max_depth=2,
                      allow_neutral=False,
                      find_all_paths=False,
                      output_file="genotypes3.csv",
                      return_output=True)
    assert os.path.exists("genotypes3.csv")
    assert issubclass(type(out3),pd.DataFrame)
    assert len(out3) < len(out2)
    
    # No output file
    os.remove("genotypes1.csv")
    os.remove("genotypes2.csv")
    os.remove("genotypes3.csv")
    gc = copy.deepcopy(ens_with_fitness_two_site["gc"])
    out4 = pathfinder(gc=gc,
                      max_depth=2,
                      allow_neutral=False,
                      find_all_paths=False,
                      output_file=None,
                      return_output=True)
    assert len(glob.glob("*.csv")) == 0
    assert issubclass(type(out4),pd.DataFrame)

    # No output, write output
    gc = copy.deepcopy(ens_with_fitness_two_site["gc"])
    out5 = pathfinder(gc=gc,
                      max_depth=2,
                      allow_neutral=False,
                      find_all_paths=False,
                      output_file="genotypes.csv",
                      return_output=False)
    assert os.path.exists("genotypes.csv")
    assert out5 is None

    # test find_all_paths
    # No output, write output
    gc = copy.deepcopy(ens_with_fitness_two_site["gc"])
    not_all = pathfinder(gc=gc,
                         max_depth=2,
                         allow_neutral=True,
                         find_all_paths=False,
                         output_file=None,
                         return_output=True)
    
    gc = copy.deepcopy(ens_with_fitness_two_site["gc"])
    all_path = pathfinder(gc=gc,
                         max_depth=2,
                         allow_neutral=True,
                         find_all_paths=True,
                         output_file=None,
                         return_output=True)

    assert len(not_all) < len(all_path)

    # Pass in bad values
    gc = copy.deepcopy(ens_with_fitness_two_site["gc"])
    for v in variable_types["everything"]:
        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            pathfinder(gc=v,max_depth=1,allow_neutral=True)

    gc = copy.deepcopy(ens_with_fitness_two_site["gc"])
    for v in variable_types["not_ints_or_coercable"]:
        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            pathfinder(gc=gc,max_depth=v,allow_neutral=True)

    for v in [-1,0]:
        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            pathfinder(gc=gc,max_depth=v,allow_neutral=True)

    gc = copy.deepcopy(ens_with_fitness_two_site["gc"])
    for v in variable_types["not_bools"]:
        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            pathfinder(gc=gc,max_depth=1,allow_neutral=v)
    
    gc = copy.deepcopy(ens_with_fitness_two_site["gc"])
    for v in variable_types["not_bools"]:
        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            pathfinder(gc=gc,max_depth=1,find_all_paths=v)

    gc = copy.deepcopy(ens_with_fitness_two_site["gc"])
    for v in variable_types["not_bools"]:
        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            pathfinder(gc=gc,max_depth=1,return_output=v)

    os.chdir(current_dir)


