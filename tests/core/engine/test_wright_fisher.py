import pytest

from eee.core.engine.wright_fisher import write_wf_outputs
from eee.core.engine.wright_fisher import wright_fisher

from eee.core.genotype import Genotype

import numpy as np
import pandas as pd

import os
import glob
import copy
import pickle

def test_write_wf_outputs(ens_with_fitness,tmpdir):
    
    def _check_write_integrity(pickle_file,generations):

        generations = [dict(zip(*g)) for g in generations]

        with open(pickle_file,'rb') as f:
            from_file = pickle.load(f)

        assert len(from_file) == len(generations)
        for i in range(len(generations)):
            for k in generations[i]:
                assert generations[i][k] == from_file[i][k]
                

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    # Run short simulation to use to write outputs
    gc = ens_with_fitness["gc"]
    template_gc, generations = wright_fisher(gc=gc,
                                             mutation_rate=0.001,
                                             population=10000,
                                             num_generations=1000)
    
    template_generations = []
    for g in generations:
        a = np.array(list(g.keys()))
        b = np.array(list(g.values()))
        template_generations.append((a,b))
    
    # --------------------------------------------------------------------------
    # Write prefix is None -- do not write out

    gc = copy.deepcopy(template_gc)
    generations = copy.deepcopy(template_generations)

    out_gc, out_gen = write_wf_outputs(gc=gc,
                                     generations=generations,
                                     write_prefix=None,
                                     write_counter=0,
                                     num_write_digits=3,
                                     final_dump=False)
    
    assert out_gc is gc
    assert out_gen is generations
    assert len(out_gc.fitnesses) == len(gc.fitnesses)
    assert len(out_gen) == len(generations)
    assert len(glob.glob("*.pickle")) == 0
    assert len(glob.glob("*.csv")) == 0

    # --------------------------------------------------------------------------
    # Write out non-final file

    gc = copy.deepcopy(template_gc)
    generations = copy.deepcopy(template_generations)

    total_num_genotypes = len(gc.genotypes)

    out_gc, out_gen = write_wf_outputs(gc=gc,
                                     generations=generations,
                                     write_prefix="test",
                                     write_counter=0,
                                     num_write_digits=3,
                                     final_dump=False)
    
    assert out_gc is gc
    assert out_gen is not generations
    assert len(out_gc.fitnesses) == len(gc.fitnesses)
    assert len(out_gen) == 1
    assert len(glob.glob("*.pickle")) == 1
    
    _check_write_integrity("test_generations_000.pickle",
                           template_generations[:-1])
        
    for f in glob.glob("*.pickle"):
        os.remove(f)

    # Make sure we wrote out the right number of genotypes
    df = pd.read_csv("test_genotypes.csv")
    assert total_num_genotypes - len(gc.genotypes) == len(df)
    in_gc = set(list(gc.genotypes.keys()))
    in_out = set(list(df.genotype))
    assert len(list(in_gc.intersection(in_out))) == 0
    os.remove("test_genotypes.csv")


    # --------------------------------------------------------------------------
    # Write out final file

    gc = copy.deepcopy(template_gc)
    generations = copy.deepcopy(template_generations)

    out_gc, out_gen = write_wf_outputs(gc=gc,
                                     generations=generations,
                                     write_prefix="test",
                                     write_counter=0,
                                     num_write_digits=3,
                                     final_dump=True)
    
    assert out_gc is gc
    assert out_gen is not generations
    assert len(out_gc.fitnesses) == len(gc.fitnesses)
    assert len(out_gen) == 0
    assert len(glob.glob("*.pickle")) == 1
    
    _check_write_integrity("test_generations_000.pickle",
                           template_generations)
        
    for f in glob.glob("*.pickle"):
        os.remove(f)

    # Make sure we wrote out the right number of genotypes
    df = pd.read_csv("test_genotypes.csv")
    assert total_num_genotypes == len(df)
    os.remove("test_genotypes.csv")

    # --------------------------------------------------------------------------
    # Write out final file, altered write_counter

    gc = copy.deepcopy(template_gc)
    generations = copy.deepcopy(template_generations)

    out_gc, out_gen = write_wf_outputs(gc=gc,
                                     generations=generations,
                                     write_prefix="test",
                                     write_counter=10,
                                     num_write_digits=3,
                                     final_dump=True)
    
    assert out_gc is gc
    assert out_gen is not generations
    assert len(out_gc.fitnesses) == len(gc.fitnesses)
    assert len(out_gen) == 0
    assert len(glob.glob("*.pickle")) == 1
    
    _check_write_integrity("test_generations_010.pickle",
                           template_generations)
    for f in glob.glob("*.pickle"):
        os.remove(f)

    # Make sure we wrote out the right number of genotypes
    df = pd.read_csv("test_genotypes.csv")
    assert total_num_genotypes == len(df)
    os.remove("test_genotypes.csv")

    # --------------------------------------------------------------------------
    # Write out final file, altered write_digits

    gc = copy.deepcopy(template_gc)
    generations = copy.deepcopy(template_generations)

    out_gc, out_gen = write_wf_outputs(gc=gc,
                                     generations=generations,
                                     write_prefix="test",
                                     write_counter=10,
                                     num_write_digits=2,
                                     final_dump=True)
    
    assert out_gc is gc
    assert out_gen is not generations
    assert len(out_gc.fitnesses) == len(gc.fitnesses)
    assert len(out_gen) == 0
    assert len(glob.glob("*.pickle")) == 1
    
    _check_write_integrity("test_generations_10.pickle",
                           template_generations)
    for f in glob.glob("*.pickle"):
        os.remove(f)

    os.chdir(current_dir)


def test_wright_fisher(ens_test_data,ens_with_fitness,variable_types,tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)


    def _count_pop(generation):
        return sum([generation[s] for s in generation])

    ens = ens_test_data["ens"]
    fitness_function = ens_test_data["fc"].fitness
    ddg_df = ens_test_data["ddg_df"]

    # --------------------------------------------------------------------------
    # Test three ways of setting the population (single number, dictionary, list)

    gc = Genotype(ens=ens,
                  fitness_function=fitness_function,
                  ddg_df=ddg_df)
    gc, generations = wright_fisher(gc,
                                    population=10,
                                    mutation_rate=0.1,
                                    num_generations=2)
    
    assert len(generations) == 2
    assert _count_pop(generations[0]) == 10
    assert _count_pop(generations[1]) == 10
    assert len(generations[0]) == 1
    assert generations[0][0] == 10
    assert len(glob.glob("*.pickle")) == 0 
    assert len(glob.glob("*.csv")) == 0 

    gc = Genotype(ens=ens,
                  fitness_function=fitness_function,
                  ddg_df=ddg_df)
    gc, generations = wright_fisher(gc,
                                    population={0:10},
                                    mutation_rate=0.1,
                                    num_generations=2)
    
    assert len(generations) == 2
    assert _count_pop(generations[0]) == 10
    assert _count_pop(generations[1]) == 10
    assert len(generations[0]) == 1
    assert generations[0][0] == 10
    assert len(glob.glob("*.pickle")) == 0 
    assert len(glob.glob("*.csv")) == 0 

    gc = Genotype(ens=ens,
                  fitness_function=fitness_function,
                  ddg_df=ddg_df)
    gc, generations = wright_fisher(gc,
                                    population=np.zeros(10,dtype=int),
                                    mutation_rate=0.1,
                                    num_generations=2)
    
    assert len(generations) == 2
    assert _count_pop(generations[0]) == 10
    assert _count_pop(generations[1]) == 10
    assert len(generations[0]) == 1
    assert generations[0][0] == 10
    assert len(glob.glob("*.pickle")) == 0 
    assert len(glob.glob("*.csv")) == 0 

    # --------------------------------------------------------------------------
    # Mutation rate

    gc = Genotype(ens=ens,
                  fitness_function=fitness_function,
                  ddg_df=ddg_df)
    gc, generations = wright_fisher(gc,
                                    population=10,
                                    mutation_rate=1,
                                    num_generations=2)
    high_mut_rate_genotypes = len(gc.genotypes)

    gc = Genotype(ens=ens,
                  fitness_function=fitness_function,
                  ddg_df=ddg_df)
    gc, generations = wright_fisher(gc,
                                    population=10,
                                    mutation_rate=0.0001,
                                    num_generations=2)
    low_mut_rate_genotypes = len(gc.genotypes)

    assert high_mut_rate_genotypes > low_mut_rate_genotypes
    assert len(glob.glob("*.pickle")) == 0 
    assert len(glob.glob("*.csv")) == 0 

    # --------------------------------------------------------------------------
    # Num generations

    gc = Genotype(ens=ens,
                  fitness_function=fitness_function,
                  ddg_df=ddg_df)
    gc, generations = wright_fisher(gc,
                                    population=10,
                                    mutation_rate=0.1,
                                    num_generations=10)
    assert len(generations) == 10
    assert len(glob.glob("*.pickle")) == 0 
    assert len(glob.glob("*.csv")) == 0 

    # pass in special case with 0 generations run. 
    gc = Genotype(ens=ens,
                  fitness_function=fitness_function,
                  ddg_df=ddg_df)
    gc, generations = wright_fisher(gc,
                                    population=10,
                                    mutation_rate=0.1,
                                    num_generations=0)
    assert len(generations) == 1
    assert len(glob.glob("*.pickle")) == 0 
    assert len(glob.glob("*.csv")) == 0 

    # --------------------------------------------------------------------------
    # num_mutations

    # 1 mutation
    gc = Genotype(ens=ens,
                  fitness_function=fitness_function,
                  ddg_df=ddg_df)
    gc, generations = wright_fisher(gc,
                                    population=100,
                                    mutation_rate=0.01,
                                    num_generations=1000,
                                    num_mutations=1)
    assert len(generations) < 1000
    to_sort = []
    keys = list(generations[-1].keys())
    for k in keys:
        to_sort.append((generations[-1][k],k))
    to_sort.sort()
    assert len(gc.genotypes[to_sort[-1][1]].mutations_accumulated) >= 1
    assert len(glob.glob("*.pickle")) == 0
    assert len(glob.glob("*.csv")) == 0  

    # 2 mutations
    gc = Genotype(ens=ens,
                  fitness_function=fitness_function,
                  ddg_df=ddg_df)
    gc, generations = wright_fisher(gc,
                                    population=100,
                                    mutation_rate=0.01,
                                    num_generations=1000,
                                    num_mutations=2)
    assert len(generations) < 1000
    to_sort = []
    keys = list(generations[-1].keys())
    for k in keys:
        to_sort.append((generations[-1][k],k))
    to_sort.sort()
    assert len(gc.genotypes[to_sort[-1][1]].mutations_accumulated) >= 2
    assert len(glob.glob("*.pickle")) == 0 
    assert len(glob.glob("*.csv")) == 0 

    # 3 mutations
    gc = Genotype(ens=ens,
                  fitness_function=fitness_function,
                  ddg_df=ddg_df)
    gc, generations = wright_fisher(gc,
                                    population=100,
                                    mutation_rate=0.01,
                                    num_generations=1000,
                                    num_mutations=3)
    assert len(generations) < 1000
    to_sort = []
    keys = list(generations[-1].keys())
    for k in keys:
        to_sort.append((generations[-1][k],k))
    to_sort.sort()
    assert len(gc.genotypes[to_sort[-1][1]].mutations_accumulated) >= 3
    assert len(glob.glob("*.pickle")) == 0 

    # Too short to complete 50 mutations -- should run 10 generations and warn
    gc = Genotype(ens=ens,
                  fitness_function=fitness_function,
                  ddg_df=ddg_df)
    with pytest.warns():
        gc, generations = wright_fisher(gc,
                                        population=10,
                                        mutation_rate=0.001,
                                        num_generations=10,
                                        num_mutations=100)
    assert len(generations) == 10
    assert len(glob.glob("*.pickle")) == 0
    assert len(glob.glob("*.csv")) == 0  

    # --------------------------------------------------------------------------
    # verbose

    # this is not a real test, but at least will test ability for the code to 
    # run without crashing when this is set. 

    gc = Genotype(ens=ens,
                  fitness_function=fitness_function,
                  ddg_df=ddg_df)
    gc, generations = wright_fisher(gc,
                                    population=10,
                                    mutation_rate=0.1,
                                    num_generations=10,
                                    verbose=False)


    # --------------------------------------------------------------------------
    # Check the sorts of mutations being generated

    gc = Genotype(ens=ens,
                  fitness_function=fitness_function,
                  ddg_df=ddg_df)
    gc, generations = wright_fisher(gc,
                                    population=10,
                                    mutation_rate=1,
                                    num_generations=2)
    high_mut_rate_genotypes = len(gc.genotypes)

    gc = Genotype(ens=ens,
                  fitness_function=fitness_function,
                  ddg_df=ddg_df)
    gc, generations = wright_fisher(gc,
                                    population=10,
                                    mutation_rate=0.0001,
                                    num_generations=2)
    low_mut_rate_genotypes = len(gc.genotypes)

    assert high_mut_rate_genotypes > low_mut_rate_genotypes
    assert len(glob.glob("*.pickle")) == 0
    assert len(glob.glob("*.csv")) == 0 

    # --------------------------------------------------------------------------
    # Check variables

    gc = Genotype(ens=ens,
                  fitness_function=fitness_function,
                  ddg_df=ddg_df)

    not_allowed = variable_types["not_ints_or_coercable"]
    for v in not_allowed:
        print(v,type(v),flush=True)

        with pytest.raises(ValueError):
            gc, generations = wright_fisher(gc,
                                            population=v,
                                            mutation_rate=0.1,
                                            num_generations=2)

    not_allowed = variable_types["everything"]
    for v in not_allowed:
        print(v,type(v),flush=True)

        with pytest.raises(ValueError):
            gc, generations = wright_fisher(gc=v,
                                            population=10,
                                            mutation_rate=0.1,
                                            num_generations=2)


    gc = Genotype(ens=ens,
                  fitness_function=fitness_function,
                  ddg_df=ddg_df)

    not_allowed = variable_types["everything"]
    for v in not_allowed:

        if v is None:
            continue

        print(v,type(v),flush=True)

        with pytest.raises(ValueError):
            gc, generations = wright_fisher(gc,
                                            population=10,
                                            mutation_rate=0.1,
                                            num_generations=2,
                                            rng=v)

    # Check reproducibility after passing in a fixed seed (0)
    rng = np.random.Generator(np.random.PCG64(0))
    gc, generations = wright_fisher(gc,
                                    population=10,
                                    mutation_rate=0.1,
                                    num_generations=5,
                                    rng=rng)

    assert generations[0][0] == 10
    assert generations[1][0] == 9
    assert generations[1][1] == 1
    assert generations[2][0] == 8
    assert generations[2][2] == 1
    assert generations[2][3] == 1
    assert generations[3][0] == 8
    assert generations[3][3] == 2
    assert generations[4][0] == 7
    assert generations[4][3] == 3

    # --------------------------------------------------------------------------
    # Run a short calculation with a dataset that has three genotypes, one of
    # which is very fit, one that is very unfit, and wildtype which is in the 
    # middle

    gc = copy.deepcopy(ens_with_fitness["gc"])

    gc, generations = wright_fisher(gc=gc,
                                    mutation_rate=0.001,
                                    population=10000,
                                    num_generations=1000)

    genotypes = {}
    for gen in generations:
        for g in gen:
            mut_genotype = "".join(gc.genotypes[g].mutations)
            if mut_genotype not in genotypes:
                genotypes[mut_genotype] = 0
            
            genotypes[mut_genotype] += gen[g]


    assert genotypes["A1P"] < 10000
    assert genotypes[""] > genotypes["A1P"]
    assert genotypes[""] < genotypes["A1V"]
    assert genotypes["A1V"] > 9000000
    assert len(glob.glob("*.pickle")) == 0 
    assert len(glob.glob("*.csv")) == 0 


    # --------------------------------------------------------------------------
    # Check ability to write outputs

    gc = copy.deepcopy(ens_with_fitness["gc"])

    gc, generations = wright_fisher(gc=gc,
                                    mutation_rate=0.001,
                                    population=10,
                                    num_generations=1000,
                                    write_prefix="test",
                                    write_frequency=100000)
    assert len(glob.glob("*.pickle")) == 1
    assert glob.glob("*.pickle")[0] == "test_generations_0.pickle"
    assert len(generations) == 0
    for f in glob.glob("*.pickle"):
        os.remove(f)

    # Change prefix
    gc = copy.deepcopy(ens_with_fitness["gc"])

    gc, generations = wright_fisher(gc=gc,
                                    mutation_rate=0.001,
                                    population=10,
                                    num_generations=1000,
                                    write_prefix="test_this",
                                    write_frequency=100000)
    assert len(glob.glob("*.pickle")) == 1
    assert glob.glob("*.pickle")[0] == "test_this_generations_0.pickle"
    assert len(generations) == 0
    for f in glob.glob("*.pickle"):
        os.remove(f)

    # Change write frequency.
    gc = copy.deepcopy(ens_with_fitness["gc"])

    gc, generations = wright_fisher(gc=gc,
                                    mutation_rate=0.001,
                                    population=10,
                                    num_generations=1000,
                                    write_prefix="test",
                                    write_frequency=100)
    assert len(glob.glob("*.pickle")) == 10
    files = list(glob.glob("*.pickle"))
    files.sort()
    assert files[0] == "test_generations_00.pickle"
    assert files[9] == "test_generations_09.pickle"
    assert len(generations) == 0
    for f in glob.glob("*.pickle"):
        os.remove(f)

    os.chdir(current_dir)
