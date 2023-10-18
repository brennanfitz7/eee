import pytest

from eee.core.engine.follow_tree import _simulate_branch
from eee.core.engine.follow_tree import follow_tree
from eee.core.genotype.single_genotype import SingleGenotype
from eee.core.genotype import Genotype

from eee.io import read_tree

import numpy as np
import pandas as pd
import ete3

import copy
import os
import pickle
import glob

def test__simulate_branch(ens_with_fitness,newick_files,tmpdir):

    def _get_fresh_nodes(newick_file):

        tree = read_tree(newick_file)

        for n in tree.traverse(strategy="levelorder"):
            if not n.is_leaf():
                start = n.copy()
                end, _ = start.get_children()
                break

        start.add_feature("population",{0:50})
        assert hasattr(start,"population")
        assert not hasattr(end,"population")

        start.name = "anc00"
        end.name = "anc01"

        return tree, start, end

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    template_gc = ens_with_fitness["gc"]

    gc = copy.deepcopy(template_gc)    
    rng = np.random.Generator(np.random.PCG64(1))
    _, start, end = _get_fresh_nodes(newick_files["simple.newick"])

    _simulate_branch(start_node=start,
                     end_node=end,
                     gc=gc,
                     mutation_rate=0.1,
                     num_generations=1000,
                     write_prefix="test",
                     rng=rng)

    assert hasattr(end,"population")
    assert issubclass(type(end.population),dict)
    for g in end.population:
        assert issubclass(type(gc.genotypes[g]),SingleGenotype)

    assert start.name == "anc00"
    assert end.name == "anc01"
    assert os.path.exists("test_anc00-anc01.pickle")
    with open("test_anc00-anc01.pickle","rb") as f:
        gens = pickle.load(f)
    assert issubclass(type(gens),list)
    assert issubclass(type(gens[0]),dict)
    assert gens[0] == start.population
    assert gens[-1] == end.population
    assert len(gens) < 1001
    assert gens[0][0] == 50
    assert gens[-1][0] < 50
    for f in glob.glob("test*.pickle"):
        os.remove(f)

    gc = copy.deepcopy(template_gc)    
    rng = np.random.Generator(np.random.PCG64(1))
    _, start, end = _get_fresh_nodes(newick_files["simple.newick"])

    # Too short to get all of the requested mutations
    with pytest.warns():
        _simulate_branch(start_node=start,
                        end_node=end,
                        gc=gc,
                        mutation_rate=0.1,
                        num_generations=2,
                        write_prefix="test",
                        rng=rng)
    assert os.path.exists("test_anc00-anc01.pickle")
    with open("test_anc00-anc01.pickle","rb") as f:
        gens = pickle.load(f)
    assert issubclass(type(gens),list)
    assert issubclass(type(gens[0]),dict)
    assert gens[0] == start.population
    assert gens[-1] == end.population
    assert len(gens) == 2
    assert gens[0][0] == 50
    assert gens[-1][0] < 50
    for f in glob.glob("test*.pickle"):
        os.remove(f)

    gc = copy.deepcopy(template_gc)    
    rng = np.random.Generator(np.random.PCG64(1))
    _, start, end = _get_fresh_nodes(newick_files["simple.newick"])

    # Too short to get all of the requested mutations. Turn mutation rate 
    # waaaaay down and make sure it manifests as no mutations occuring in a 
    # single generation. 
    with pytest.warns():
        _simulate_branch(start_node=start,
                        end_node=end,
                        gc=gc,
                        mutation_rate=1e-100,
                        num_generations=2,
                        write_prefix="test",
                        rng=rng)
    assert os.path.exists("test_anc00-anc01.pickle")
    with open("test_anc00-anc01.pickle","rb") as f:
        gens = pickle.load(f)
    assert issubclass(type(gens),list)
    assert issubclass(type(gens[0]),dict)
    assert gens[0] == start.population
    assert gens[-1] == end.population
    assert len(gens) == 2
    assert gens[0][0] == 50
    assert gens[-1][0] == 50
    for f in glob.glob("test*.pickle"):
        os.remove(f)

    os.chdir(current_dir)


def test_follow_tree(ens_with_fitness,newick_files,variable_types,tmpdir):
    
    current_dir = os.getcwd()
    os.chdir(tmpdir)

    template_gc = ens_with_fitness["gc"]
    tree = ete3.Tree(newick_files["simple.newick"])
    rng = np.random.Generator(np.random.PCG64(1))

    gc = copy.deepcopy(template_gc)  

    assert len(glob.glob("*.pickle")) == 0
    gc_out, tree_out = follow_tree(gc,
                                   tree=tree,
                                   population=100,
                                   mutation_rate=0.1,
                                   num_generations=1000,
                                   burn_in_generations=10,
                                   write_prefix="eee_tree",
                                   rng=rng)
    assert issubclass(type(gc_out),Genotype)
    assert issubclass(type(tree_out),ete3.TreeNode)
    
    all_genotypes_seen = []
    expected_files = ["eee_tree.newick",
                      "eee_tree_anc01-A.pickle",
                      "eee_tree_anc02-D.pickle",
                      "eee_tree_anc00-anc01.pickle",
                      "eee_tree_anc01-B.pickle",
                      "eee_tree_burn-in-anc00.pickle",
                      "eee_tree_anc00-anc02.pickle",
                      "eee_tree_anc02-C.pickle",
                      "eee_tree_genotypes.csv"]
    for f in expected_files:
        assert os.path.isfile(f)

        if f.split(".")[-1] == "pickle":

            with open(f,'rb') as handle:
                gens = pickle.load(handle)
            assert issubclass(type(gens),list)
            assert len(gens) >= 2
            assert issubclass(type(gens[0]),dict)

            for g in gens:
                all_genotypes_seen.extend(g.keys())

        if f.split(".")[-1] == "csv":
            out_df = pd.read_csv(f) 
            # just make sure has right columns
            assert "genotype" in out_df.columns
            assert "trajectory" in out_df.columns

        if f.split(".")[-1] == "newick":
            # make sure we can read it
            out_tree = read_tree(f)
            
    all_genotypes_seen = set(all_genotypes_seen)
    in_df = set(out_df["genotype"])
    assert all_genotypes_seen == in_df
    
    # Make sure generations are passing one to the other properly
    def _read_and_compare(pickle1,pickle2):

        with open(pickle1,'rb') as f:
            start = pickle.load(f)
    
        with open(pickle2,'rb') as f:
            end = pickle.load(f)

        assert start[-1] == end[0]

        # make sure check is actually happening
        assert start[0] != end[-1]

    #
    _read_and_compare("eee_tree_burn-in-anc00.pickle",
                      "eee_tree_anc00-anc01.pickle")

    _read_and_compare("eee_tree_burn-in-anc00.pickle",
                      "eee_tree_anc00-anc02.pickle")

    _read_and_compare("eee_tree_anc00-anc01.pickle",
                      "eee_tree_anc01-A.pickle")

    _read_and_compare("eee_tree_anc00-anc01.pickle",
                      "eee_tree_anc01-B.pickle")  

    _read_and_compare("eee_tree_anc00-anc02.pickle",
                      "eee_tree_anc02-C.pickle")

    _read_and_compare("eee_tree_anc00-anc02.pickle",
                      "eee_tree_anc02-D.pickle")  
    
    # test the test
    with pytest.raises(AssertionError):
        _read_and_compare("eee_tree_burn-in-anc00.pickle",
                          "eee_tree_anc02-D.pickle")  

    for f in glob.glob("*.*"):
        os.remove(f)

    # Make sure run happens with default args before passing in bad ones
    tree = ete3.Tree(newick_files["simple.newick"])
    gc = copy.deepcopy(template_gc)
    follow_tree(gc=gc,
                num_generations=1000,
                mutation_rate=0.1,
                tree=tree)
    for f in glob.glob("*.*"):
        os.remove(f)

    tree = ete3.Tree(newick_files["simple.newick"])
    gc = copy.deepcopy(template_gc) 
    for v in variable_types["everything"]:
        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            follow_tree(gc=v,
                        num_generations=1000,
                        mutation_rate=0.1,
                        tree=tree)
    for f in glob.glob("*.*"):
        os.remove(f)

    tree = ete3.Tree(newick_files["simple.newick"])
    gc = copy.deepcopy(template_gc) 
    for v in variable_types["everything"]:
        if v is None:
            continue
        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            follow_tree(gc=gc,
                        num_generations=1000,
                        mutation_rate=0.1,
                        tree=v)
    for f in glob.glob("*.*"):
        os.remove(f)

    tree = ete3.Tree(newick_files["simple.newick"])
    gc = copy.deepcopy(template_gc) 
    for v in variable_types["not_ints_or_coercable"]:
        if hasattr(v,"__iter__"):

            bad_types = [pd.DataFrame,type]
            skip = True
            for b in bad_types:
                if issubclass(type(v),b):
                    skip = False
            
            if skip:
                continue

        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            follow_tree(gc=gc,
                        population=v,
                        num_generations=1000,
                        mutation_rate=0.1,
                        tree=tree)
    for f in glob.glob("*.*"):
        os.remove(f)

    tree = ete3.Tree(newick_files["simple.newick"])
    gc = copy.deepcopy(template_gc) 
    for v in [-1,0]:
        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            follow_tree(gc=gc,
                        population=v,
                        num_generations=1000,
                        mutation_rate=0.1,
                        tree=tree)
    for f in glob.glob("*.*"):
        os.remove(f)

    tree = ete3.Tree(newick_files["simple.newick"])
    gc = copy.deepcopy(template_gc) 
    for v in variable_types["not_ints_or_coercable"]:
        
        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            follow_tree(gc=gc,
                        num_generations=v,
                        mutation_rate=0.1,
                        tree=tree)
    for f in glob.glob("*.*"):
        os.remove(f)

    tree = ete3.Tree(newick_files["simple.newick"])
    gc = copy.deepcopy(template_gc) 
    for v in [-1]:
        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            follow_tree(gc=gc,
                        num_generations=v,
                        mutation_rate=0.1,
                        tree=tree)
    for f in glob.glob("*.*"):
        os.remove(f)

    tree = ete3.Tree(newick_files["simple.newick"])
    gc = copy.deepcopy(template_gc) 
    for v in variable_types["not_floats_or_coercable"]:
        
        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            follow_tree(gc=gc,
                        num_generations=1000,
                        mutation_rate=v,
                        tree=tree)
    for f in glob.glob("*.*"):
        os.remove(f)

    tree = ete3.Tree(newick_files["simple.newick"])
    gc = copy.deepcopy(template_gc) 
    for v in [-1,0]:
        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            follow_tree(gc=gc,
                        num_generations=1000,
                        mutation_rate=v,
                        tree=tree)
    for f in glob.glob("*.*"):
        os.remove(f)

    tree = ete3.Tree(newick_files["simple.newick"])
    gc = copy.deepcopy(template_gc) 
    for v in variable_types["not_ints_or_coercable"]:
        
        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            follow_tree(gc=gc,
                        num_generations=1000,
                        mutation_rate=0.1,
                        burn_in_generations=v,
                        tree=tree)
    for f in glob.glob("*.*"):
        os.remove(f)

    tree = ete3.Tree(newick_files["simple.newick"])
    gc = copy.deepcopy(template_gc) 
    for v in [-1]:
        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            follow_tree(gc=gc,
                        num_generations=1000,
                        mutation_rate=0.1,
                        burn_in_generations=v,
                        tree=tree)
    for f in glob.glob("*.*"):
        os.remove(f)

    tree = ete3.Tree(newick_files["simple.newick"])
    gc = copy.deepcopy(template_gc) 
    for v in variable_types["everything"]:
        if v is None:
            continue
        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            follow_tree(gc=gc,
                        num_generations=1000,
                        mutation_rate=0.1,
                        tree=tree,
                        rng=v)
    for f in glob.glob("*.*"):
        os.remove(f)

    tree = ete3.Tree(newick_files["simple.newick"])
    gc = copy.deepcopy(template_gc)
    with pytest.raises(ValueError):
        follow_tree(gc=gc,
                num_generations=1000,
                mutation_rate=0.1,
                write_prefix=None,
                tree=tree)

    os.chdir(current_dir)