
import pytest

from eee.core.genotype.single_genotype import SingleGenotype
from eee.core.genotype.genotype import Genotype

import numpy as np
import pandas as pd

import os
import random
import copy

def test_Genotype(ens_test_data):

    ens = ens_test_data["ens"]
    fitness_function = ens_test_data["fc"].fitness
    ddg_dict = ens_test_data["ddg_dict"]
    ddg_df = ens_test_data["ddg_df"]
    choice_function = random.choice

    gc = Genotype(ens=ens,
                           fitness_function=fitness_function,
                           ddg_df=ddg_df,
                           choice_function=choice_function)

    # Make sure that ddg_df loaded correctly into ddg_dict
    for a in ddg_dict:
        for b in ddg_dict[a]:
            assert np.array_equal(gc._ddg_dict[a][b], ddg_dict[a][b])
            
    # Make sure possible sites and mutations at sites are correct
    assert np.array_equal(gc._possible_sites,[1,2])
    assert np.array_equal(gc._mutations_at_sites[1],["M1A","M1V"])
    assert np.array_equal(gc._mutations_at_sites[2],["P2R","P2Q"])

    # Make sure it created correct SingleGenotype instance
    assert len(gc.genotypes) == 1
    assert issubclass(type(gc.genotypes[0]),SingleGenotype)
    assert gc.genotypes[0]._ens is ens
    assert gc._fitness_function is fitness_function
    assert gc.genotypes[0]._ddg_dict is gc._ddg_dict
    assert len(gc.genotypes[0].sites) == 0
    assert len(gc.genotypes[0].mutations) == 0
    assert len(gc.genotypes[0].mutations_accumulated) == 0
    assert len(gc.genotypes[0].mut_energy) == 2
    assert gc.genotypes[0].mut_energy[0] == 0
    assert gc.genotypes[0].mut_energy[1] == 0

    # Make sure trajectories and fitnesses are set up correctly. 
    assert len(gc.trajectories) == 1
    assert len(gc.trajectories[0]) == 1
    assert gc.trajectories[0][0] == 0
    assert len(gc.fitnesses) == 1
    assert gc.fitnesses[0] >= 0
    assert gc.fitnesses[0] <= 1

    assert gc._choice_function is random.choice
    gc = Genotype(ens=ens,
                           fitness_function=fitness_function,
                           ddg_df=ddg_df,
                           choice_function=None)
    assert gc._choice_function is np.random.choice


    # send in copy of ddg_df with mangled column names. Should throw a 
    # ValueError  
    bad_ddg_df = ddg_df.copy()
    new_columns = list(bad_ddg_df.columns)
    idx = new_columns.index("s1")
    new_columns[idx] = "not_a_species"
    bad_ddg_df.columns = new_columns

    with pytest.raises(ValueError):
        gc = Genotype(ens=ens,
                               fitness_function=fitness_function,
                               ddg_df=bad_ddg_df)

def test_Genotype__create_ddg_dict(ens_test_data):

    ens = ens_test_data["ens"]
    fitness_function = ens_test_data["fc"].fitness
    ddg_dict = ens_test_data["ddg_dict"]
    ddg_df = ens_test_data["ddg_df"]

    gc = Genotype(ens=ens,
                           fitness_function=fitness_function,
                           ddg_df=ddg_df)
    
    # Make sure that ddg_df loaded correctly into ddg_dict
    for a in ddg_dict:
        for b in ddg_dict[a]:
            assert np.array_equal(gc._ddg_dict[a][b], ddg_dict[a][b])

    # send in copy of ddg_df with mangled column names. Should throw a 
    # ValueError  
    bad_ddg_df = ddg_df.copy()
    new_columns = list(bad_ddg_df.columns)
    idx = new_columns.index("s1")
    new_columns[idx] = "not_a_species"
    bad_ddg_df.columns = new_columns

    with pytest.raises(ValueError):
        gc = Genotype(ens=ens,
                      fitness_function=fitness_function,
                      ddg_df=bad_ddg_df)


def test_Genotype__add_genotype(ens_test_data):
    ## Tested implicitly (and more extensively) in test_Genotype_mutate

    ens = ens_test_data["ens"]
    fitness_function = ens_test_data["fc"].fitness
    ddg_df = ens_test_data["ddg_df"]

    gc = Genotype(ens=ens,
                  fitness_function=fitness_function,
                  ddg_df=ddg_df)
    
    sc = SingleGenotype(ens=ens,
                        ddg_dict=gc.ddg_dict)
    gc._add_genotype(new_genotype=sc,
                     prev_index=0)
    
    assert gc._genotypes[1] is sc
    assert np.array_equal(gc._trajectories[1],[0,1])
    assert np.array_equal(gc._mut_energies[1],[0,0])
    assert gc._fitnesses[0] == gc._fitnesses[1]


def test_Genotype_mutate(ens_test_data):

    ens = ens_test_data["ens"]
    fitness_function = ens_test_data["fc"].fitness
    ddg_df = ens_test_data["ddg_df"]

    gc = Genotype(ens=ens,
                  fitness_function=fitness_function,
                  ddg_df=ddg_df)
    
    # Mutate residue not there. 
    with pytest.raises(IndexError):
        gc.mutate(1)

    # Mutate residue not there. 
    with pytest.raises(IndexError):
        gc.mutate(-1)

    # Chain of mutations...
    allowable_muts =set(ddg_df["mut"])
    for i in range(4):

        gc.mutate(i)
    
        # Make sure trajectories accumulating as we expect
        assert len(gc.trajectories) == i + 2
        assert len(gc.trajectories[i+1]) == i + 2
        assert np.array_equal(gc.trajectories[i+1],
                              np.arange(i+2,dtype=int))
        
        # make sure fitness values are in proper bounds
        assert len(gc.fitnesses) == i + 2
        for f in gc.fitnesses:
            assert gc.fitnesses[f] >= 0 and gc.fitnesses[f] <= 1

        # Make sure we are accumulating genotypes objects that are not the same
        assert len(gc.genotypes) == i + 2
        for j in range(len(gc.genotypes)):
            assert issubclass(type(gc.genotypes[j]),SingleGenotype)
            assert set(gc.genotypes[j].mutations).issubset(allowable_muts)

            for k in range(j+1,len(gc.genotypes)):
                assert gc.genotypes[j] is not gc.genotypes[k]

    ens = ens_test_data["ens"]
    fitness_function = ens_test_data["fc"].fitness
    ddg_df = ens_test_data["ddg_df"]

    gc = Genotype(ens=ens,
                  fitness_function=fitness_function,
                  ddg_df=ddg_df)
    
    seen_times = {"M1A":0,"M1V":0}
    for i in range(10):
        gc.mutate(index=0,site=1)
        assert gc.genotypes[i+1].sites[0] == 1
        assert len(gc.genotypes[i+1].sites) == 1

        assert len(gc.genotypes[i+1].mutations) == 1
        seen_times[gc.genotypes[i+1].mutations[0]] += 1

    assert seen_times["M1A"] > 0
    assert seen_times["M1V"] > 0

    ens = ens_test_data["ens"]
    fitness_function = ens_test_data["fc"].fitness
    ddg_df = ens_test_data["ddg_df"]

    gc = Genotype(ens=ens,
                  fitness_function=fitness_function,
                  ddg_df=ddg_df)
    
    seen_times = {"M1A":0,"M1V":0}
    for i in range(10):
        gc.mutate(index=0,site=1,mutation="M1A")
        assert gc.genotypes[i+1].sites[0] == 1
        assert len(gc.genotypes[i+1].sites) == 1

        assert len(gc.genotypes[i+1].mutations) == 1
        seen_times[gc.genotypes[i+1].mutations[0]] += 1

    assert seen_times["M1A"] == 10
    assert seen_times["M1V"] == 0

    gc = Genotype(ens=ens,
                  fitness_function=fitness_function,
                  ddg_df=ddg_df)
    
    gc.mutate(index=0,site=2)
    assert gc.genotypes[1].sites[0] == 2
    assert len(gc.genotypes[1].sites) == 1

    gc = Genotype(ens=ens,
                  fitness_function=fitness_function,
                  ddg_df=ddg_df)
    
    gc.mutate(index=0,site=2,mutation="P2R")
    assert gc.genotypes[1].sites[0] == 2
    assert len(gc.genotypes[1].sites) == 1
    assert len(gc.genotypes[1].mutations) == 1
    assert gc.genotypes[1].mutations[0] == "P2R"


def test_Genotype_conditional_mutate(ens_with_fitness):

    ens = copy.deepcopy(ens_with_fitness["ens"])
    fitness_function = copy.deepcopy(ens_with_fitness["fc"].fitness)
    
    to_df = {"site":[1,1,1],
             "mut":["A1V","A1P","A1C"],
             "s1":[-1.677,0.167,0.000],
             "s2":[ 3.333,-5000,0.000]}
    
    ddg_df = pd.DataFrame(to_df)

    gc = Genotype(ens=ens,
                  fitness_function=fitness_function,
                  ddg_df=ddg_df)
    
    # ----------------- Favorable. ----------------------------- 
    result = gc.conditional_mutate(index=0,
                                   site=1,
                                   mutation="A1V",
                                   condition_fcn=np.greater_equal)
    assert result == 1
    assert gc._genotypes[1].mutations[0] == "A1V"

    result = gc.conditional_mutate(index=0,
                                   site=1,
                                   mutation="A1V",
                                   condition_fcn=np.greater)
    assert result == 2
    assert gc._genotypes[2].mutations[0] == "A1V"

    # Not same -- do not incorporate
    result = gc.conditional_mutate(index=0,
                                   site=1,
                                   mutation="A1V",
                                   condition_fcn=np.equal)
    assert result == -1
    assert len(gc._genotypes) == 3

    # ----------------- Favorable. -----------------------------
    result = gc.conditional_mutate(index=0,
                                   site=1,
                                   mutation="A1P",
                                   condition_fcn=np.greater_equal)
    assert result == -1
    assert len(gc._genotypes) == 3

    result = gc.conditional_mutate(index=0,
                                   site=1,
                                   mutation="A1P",
                                   condition_fcn=np.greater)
    assert result == -1
    assert len(gc._genotypes) == 3

    # Not same -- do not incorporate
    result = gc.conditional_mutate(index=0,
                                   site=1,
                                   mutation="A1P",
                                   condition_fcn=np.equal)
    assert result == -1
    assert len(gc._genotypes) == 3

    # Perversely select something lowering fitness
    result = gc.conditional_mutate(index=0,
                                   site=1,
                                   mutation="A1P",
                                   condition_fcn=np.less)
    assert result == 3
    assert gc._genotypes[3].mutations[0] == "A1P"

    # ----------------- Neutral. -----------------------------
    result = gc.conditional_mutate(index=0,
                                   site=1,
                                   mutation="A1C",
                                   condition_fcn=np.greater_equal)
    assert result == 4
    assert gc._genotypes[4].mutations[0] == "A1C"

    result = gc.conditional_mutate(index=0,
                                   site=1,
                                   mutation="A1C",
                                   condition_fcn=np.greater)
    assert result == -1
    assert len(gc._genotypes) == 5

    # Not same -- do not incorporate
    result = gc.conditional_mutate(index=0,
                                   site=1,
                                   mutation="A1C",
                                   condition_fcn=np.equal)
    assert result == 5
    assert gc._genotypes[5].mutations[0] == "A1C"


def test_Genotype_dump_to_csv(ens_test_data,tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    # Generic dump
    ens = ens_test_data["ens"]
    fitness_function = ens_test_data["fc"].fitness
    ddg_df = ens_test_data["ddg_df"]

    gc = Genotype(ens=ens,
                           fitness_function=fitness_function,
                           ddg_df=ddg_df)
    
    gc.dump_to_csv(filename="test.csv",
                   keep_genotypes=None)

    df = pd.read_csv("test.csv")
    assert len(df) == 1
    os.remove("test.csv")

    assert len(gc.genotypes) == 0
    assert len(gc.trajectories) == 0
    assert len(gc.mut_energies) == 0
    assert len(gc.fitnesses) == 0

    # Dump after making 10 mutations

    gc = Genotype(ens=ens,
                           fitness_function=fitness_function,
                           ddg_df=ddg_df)
    for i in range(10):
        gc.mutate(i)
    
    gc.dump_to_csv(filename="test.csv",
                   keep_genotypes=None)

    df = pd.read_csv("test.csv")
    assert len(df) == 11
    os.remove("test.csv")

    assert len(gc.genotypes) == 0
    assert len(gc.trajectories) == 0
    assert len(gc.mut_energies) == 0
    assert len(gc.fitnesses) == 0

    # Dump only half of the genotypes

    gc = Genotype(ens=ens,
                           fitness_function=fitness_function,
                           ddg_df=ddg_df)
    for i in range(10):
        gc.mutate(i)

    keep_genotypes = [0,2,4,6,8,10]
    gc.dump_to_csv(filename="test.csv",
                   keep_genotypes=keep_genotypes)

    df = pd.read_csv("test.csv")
    assert len(df) == 5
    os.remove("test.csv")

    keys = list(gc.genotypes.keys())
    keys.sort()
    assert np.array_equal(keys,keep_genotypes)

    # Dump in two steps, appending

    gc = Genotype(ens=ens,
                           fitness_function=fitness_function,
                           ddg_df=ddg_df)
    for i in range(10):
        gc.mutate(i)

    keep_genotypes = [0,2,4,6,8,10]
    gc.dump_to_csv(filename="test.csv",
                   keep_genotypes=keep_genotypes)

    df = pd.read_csv("test.csv")
    assert len(df) == 5

    keys = list(gc.genotypes.keys())
    keys.sort()
    assert np.array_equal(keys,keep_genotypes)

    gc.dump_to_csv(filename="test.csv")
    df = pd.read_csv("test.csv")
    assert len(df) == 11
    assert len(np.unique(df.genotype)) == 11
    os.remove("test.csv")

    assert len(gc.genotypes) == 0
    assert len(gc.trajectories) == 0
    assert len(gc.mut_energies) == 0
    assert len(gc.fitnesses) == 0

    os.chdir(current_dir)

def test_Genotype_to_dict(ens_test_data):
    
    ens = ens_test_data["ens"]
    fitness_function = ens_test_data["fc"].fitness
    ddg_df = ens_test_data["ddg_df"]

    gc = Genotype(ens=ens,
                  fitness_function=fitness_function,
                  ddg_df=ddg_df)
    
    out = gc.to_dict()
    assert issubclass(type(out),dict)
    assert len(out) == 1
    assert issubclass(type(out["ddg_df"]),pd.DataFrame)

def test_Genotype_df(ens_test_data):

    # not a great test of all features of dataframe. Makes sure table has right
    # columns and that mutational parent tracking works. 

    ens = ens_test_data["ens"]
    fitness_function = ens_test_data["fc"].fitness
    ddg_df = ens_test_data["ddg_df"]

    gc = Genotype(ens=ens,
                           fitness_function=fitness_function,
                           ddg_df=ddg_df)
    
    assert issubclass(type(gc.df),pd.DataFrame)

    assert np.array_equal(gc.df.columns,
                          ["genotype",
                           "mutations","num_mutations",
                           "accum_mut","num_accum_mut",
                           "parent","trajectory",
                           "s1_ddg","s2_ddg",
                           "fitness"])
    
    assert len(gc.df) == len(gc.genotypes)

    gc.mutate(0)
    gc.mutate(0)
    gc.mutate(1)
    gc.mutate(3)

    # Drop first row from comparison because null causes problems
    assert np.array_equal(list(gc.df["parent"])[1:],[0,0,1,3])
    assert gc.df.loc[3,"num_accum_mut"] == 2


def test_Genotype_wt_sequence(ens_test_data):
    
    ens = ens_test_data["ens"]
    fitness_function = ens_test_data["fc"].fitness
    ddg_df = ens_test_data["ddg_df"]

    gc = Genotype(ens=ens,
                           fitness_function=fitness_function,
                           ddg_df=ddg_df)
    
    assert gc.wt_sequence == "MP"

    ddg_df = pd.DataFrame({"site":[1,2,3],
                           "mut":["M1A","P2Q","L3A"],
                           "s1":[0,0,0],
                           "s2":[0,0,0]})
    
    
    gc = Genotype(ens=ens,
                           fitness_function=fitness_function,
                           ddg_df=ddg_df)

    assert gc.wt_sequence == "MPL"
    

def test_Genotype_genotypes(ens_test_data):
    
    ens = ens_test_data["ens"]
    fitness_function = ens_test_data["fc"].fitness
    ddg_df = ens_test_data["ddg_df"]

    gc = Genotype(ens=ens,
                           fitness_function=fitness_function,
                           ddg_df=ddg_df)
    
    gc.mutate(0)

    assert len(gc.genotypes) == 2
    assert issubclass(type(gc.genotypes[0]),SingleGenotype)
    assert issubclass(type(gc.genotypes[1]),SingleGenotype)
    assert gc.genotypes[0] is not gc.genotypes[1]

def test_Genotype_trajectories(ens_test_data):

    ens = ens_test_data["ens"]
    fitness_function = ens_test_data["fc"].fitness
    ddg_df = ens_test_data["ddg_df"]

    gc = Genotype(ens=ens,
                           fitness_function=fitness_function,
                           ddg_df=ddg_df)
    gc.mutate(0)
    assert len(gc.trajectories) == 2
    assert len(gc.trajectories[0]) == 1
    assert gc.trajectories[0][0] == 0

    assert len(gc.trajectories[1]) == 2
    assert gc.trajectories[1][0] == 0
    assert gc.trajectories[1][1] == 1

def test_Genotype_mut_energies(ens_test_data):

    ens = ens_test_data["ens"]
    fitness_function = ens_test_data["fc"].fitness
    ddg_df = ens_test_data["ddg_df"]

    gc = Genotype(ens=ens,
                           fitness_function=fitness_function,
                           ddg_df=ddg_df)
    
    gc.mutate(0)
    assert len(gc.mut_energies) == 2
    assert len(gc.mut_energies[0]) == 2
    assert np.array_equal(gc.mut_energies[0],[0,0])
    assert len(gc.mut_energies[1]) == 2
    

def test_Genotype_fitnesses(ens_test_data):
    
    ens = ens_test_data["ens"]
    fitness_function = ens_test_data["fc"].fitness
    ddg_df = ens_test_data["ddg_df"]

    gc = Genotype(ens=ens,
                           fitness_function=fitness_function,
                           ddg_df=ddg_df)
    gc.mutate(0)
    assert len(gc.fitnesses) == 2

def test_Genotype_ddg_dict(ens_test_data):

    ens = ens_test_data["ens"]
    fitness_function = ens_test_data["fc"].fitness
    ddg_df = ens_test_data["ddg_df"]
    ddg_dict = ens_test_data["ddg_dict"]

    gc = Genotype(ens=ens,
                  fitness_function=fitness_function,
                  ddg_df=ddg_df)
        
    # Make sure property works as expected. 
    for a in ddg_dict:
        for b in ddg_dict[a]:
            assert np.array_equal(gc.ddg_dict[a][b], ddg_dict[a][b])