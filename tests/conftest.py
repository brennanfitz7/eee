import pytest

# Make sure we use the non-interactive pdf backend throughout testing
import matplotlib
matplotlib.use('pdf')

from eee.core.ensemble import Ensemble
from eee.core.fitness.fitness import Fitness
from eee.core.fitness.ff import ff_on
from eee.core.fitness.ff import ff_off
from eee.core.genotype import Genotype

import numpy as np
import pandas as pd

import os
import glob

def _file_globber(*args):
    """
    Do a glob query constructed from *args relative to this directory. Return 
    a dictionary keying the basename of each file to its absolute path. 

    Example: 
    this/complicated/path/ has the files "one.csv", "two.csv", and "blah.txt".
    If args = ["this","complicated","path","*.csv"], the function will return
    {"one.csv":abspath("this/complicated/path/one.csv"),
     "two.csv":abspath("this/complicated/path/two.csv")}.
    """ 
    
    base_dir = os.path.dirname(os.path.realpath(__file__))
    search_string = os.path.join(base_dir,*args)

    file_dict = {}
    for g in glob.glob(search_string):
        key = os.path.basename(g)
        file_dict[key] = g 

    return file_dict

@pytest.fixture(scope="module")
def test_cifs():
    """
    Dictionary holding cif files for testing. 
    """

    return _file_globber("data_for_tests","test_structures","*.cif")

@pytest.fixture(scope="module")
def test_pdbs():
    """
    Dictionary holding pdb files for testing. 
    """

    return _file_globber("data_for_tests","test_structures","*.pdb")

@pytest.fixture(scope="module")
def structure_ensembles():

    base_dir = os.path.dirname(os.path.realpath(__file__))
    search_string = os.path.join(base_dir,"data_for_tests","structure_ensembles","*")

    file_dict = {}
    for g in glob.glob(search_string):
        key = os.path.basename(g)
        rcsb_files = glob.glob(os.path.join(g,"*.*"))

        file_dict[key] = rcsb_files

    return file_dict

@pytest.fixture(scope="module")
def test_ddg():

    return _file_globber("data_for_tests","test_ddg","*.csv")

@pytest.fixture(scope="module")
def spreadsheets():

    return _file_globber("data_for_tests","spreadsheets","*")

@pytest.fixture(scope="module")
def sim_json():

    return _file_globber("data_for_tests","sim_json","*.json")

@pytest.fixture(scope="module")
def conditions_input():

    return _file_globber("data_for_tests","conditions_input","*.*")

@pytest.fixture(scope="module")
def ens_test_data():

    # Basic ensemble
    ens = Ensemble(gas_constant=1)
    ens.add_species(name="s1",
                    observable=True,
                    X=1)
    ens.add_species(name="s2",
                    observable=False,
                    Y=1)
    
    # Selection conditions
    conditions = {"X":[0,1],
                  "Y":[1,0],
                  "fitness_fcn":["on","off"],
                  "select_on":["fx_obs","fx_obs"],
                  "select_on_folded":[True,True],
                  "temperature":[298.15,298.15]}

    ddg_df = pd.DataFrame({"site":[1,1,2,2],
                           "mut":["M1A","M1V","P2R","P2Q"],
                           "s1":[1,-1,0,0],
                           "s2":[-1,1,1,0]})

    # basic ddg dictionary
    ddg_dict = {}
    ddg_dict[1] = {"M1A":np.array([1,-1]),
                   "M1V":np.array([-1,1])}
    ddg_dict[2] = {"P2R":np.array([0,1]),
                   "P2Q":np.array([0,0])}

    fc = Fitness(ens=ens,
                 conditions=conditions)
    
    out = {"ens":ens,
           "conditions":conditions,
           "ddg_df":ddg_df,
           "ddg_dict":ddg_dict,
           "fc":fc}
    
    return out


@pytest.fixture(scope="module")
def ens_with_fitness():
    

    ens = Ensemble(gas_constant=1)
    ens.add_species(name="s1",
                    dG0=0,
                    observable=True)
    ens.add_species(name="s2",
                    dG0=0.167,
                    observable=False,
                    X=2)

    # Selection conditions
    conditions = {"X":[0,3.333],
                  "fitness_fcn":["on","off"],
                  "select_on":"fx_obs",
                  "select_on_folded":True,
                  "temperature":1}


    ddg_df = pd.DataFrame({"site":[1,1],
                           "mut":["A1V","A1P"],
                           "s1":[-1.677,0.167],
                           "s2":[3.333,-5000]})

    fc = Fitness(ens=ens,
                 conditions=conditions)

    gc = Genotype(ens=ens,
                  fitness_function=fc.fitness,
                  ddg_df=ddg_df)
    

    out = {"ens":ens,
           "conditions":conditions,
           "ddg_df":ddg_df,
           "fc":fc,
           "gc":gc}
    
    return out

@pytest.fixture(scope="module")
def ens_with_fitness_two_site():
    

    ens = Ensemble(gas_constant=1)
    ens.add_species(name="s1",
                    dG0=0,
                    observable=True)
    ens.add_species(name="s2",
                    dG0=0.167,
                    observable=False,
                    X=2)

    # Selection conditions
    conditions = {"X":[0,3.333],
                  "fitness_fcn":["on","off"],
                  "select_on":"fx_obs",
                  "select_on_folded":True,
                  "temperature":1}

    ddg_df = pd.DataFrame({"site":[1,1,2],
                           "mut":["A1V","A1P","A2C"],
                           "s1":[-1.677,0.167,0],
                           "s2":[3.333,-5000,0]})

    fc = Fitness(ens=ens,
                 conditions=conditions)

    gc = Genotype(ens=ens,
                  fitness_function=fc.fitness,
                  ddg_df=ddg_df)
    

    out = {"ens":ens,
           "conditions":conditions,
           "ddg_df":ddg_df,
           "fc":fc,
           "gc":gc}
    
    return out

@pytest.fixture(scope="module")
def newick_files():

    return _file_globber("data_for_tests","newick_files","*.newick")

@pytest.fixture(scope="module")
def ensemble_inputs():

    return _file_globber("data_for_tests","ensemble_inputs","*.*")

@pytest.fixture(scope="module")
def programs():
    """
    Dictionary holding paths pointing to programs to run.
    """

    dir = os.path.dirname(os.path.realpath(__file__))
    base_dir = os.path.abspath(os.path.join(dir,"data_for_tests","programs"))
    files = os.listdir(base_dir)

    out_dict = {}
    for f in files:
        out_dict[f] = os.path.join(base_dir,f)

    return out_dict

@pytest.fixture(scope="module")
def tiny_sim_output():
    """
    Tiny but complete simulation outputs.
    """

    return _file_globber("data_for_tests","tiny_sim_output","*")


@pytest.fixture(scope="module")
def variable_types():
    """
    Returns a dictionary with a bunch of different argument types to jam into
    python functions for testing. 
    """

    floats = [-1.0,0.0,1.0]
    for f in [np.float16,np.float32,np.float64]:
        for i in [-1,0,1]:
            floats.append(f(i))

    ints = [-1,0,1]
    for f in [np.int8,np.int16,np.int32,np.int64,np.intc]:
        for i in [-1,0,1]:
            ints.append(f(i))

    trues = [True,1.0,1,np.bool_(1)]
    falses = [False,0.0,0,np.bool_(0)]
    class TestClass: pass
    types = [dict,float,int,bool,str,np.ndarray,set,tuple,pd.DataFrame,TestClass]
    strings = ["test","a","",np.array(["test"])[0]]
    nones = [None]
    nans = [np.nan,pd.NA]
    infs = [-np.inf,np.inf]
    iterables = [[],(),"",{},
                 [1,2,3],(1,2,3),"abc",{1:1,2:2,3:3},
                 pd.DataFrame({"test":[1,2,3]}),
                 np.arange(3)]
    coercable_to_float = [-1,0,1,"-1","0","1","-1.0","0.0","-1.0"]
    coercable_to_int = [-1.0,0.0,1.0,"-1","0","1"]

    floats_set = set(floats)
    ints_set = set(ints)
    trues_set = set(trues)
    falses_set = set(falses)
    coercable_to_float_set = set(coercable_to_float)
    coercable_to_int_set = set(coercable_to_int)
    types_set = set(types)
    strings_set = set(strings)
    nones_set = set(nones)
    nans_set = set(nans)
    infs_set = set(infs)

    everything_set = floats_set | ints_set
    everything_set = everything_set | trues_set | falses_set
    everything_set = everything_set | coercable_to_float_set | coercable_to_int_set
    everything_set = everything_set | types_set | strings_set
    everything_set = everything_set | nones_set | nans_set | infs_set

    out = {}
    
    out["floats"] = floats
    out["ints"] = ints
    out["trues"] = trues
    out["falses"] = falses
    out["bools"] = list(trues_set | falses_set)
    out["types"] = types
    out["strings"] = strings
    out["nones"] = nones
    out["nans"] = nans
    out["infs"] = infs
    out["iterables"] = iterables
    out["coercable_to_float"] = coercable_to_float
    out["coercable_to_int"] = coercable_to_int
    out["floats_or_coercable"] = list(floats_set | coercable_to_float_set)
    out["ints_or_coercable"] = list(ints_set | coercable_to_int_set)

    out["not_floats"] = list(everything_set - floats_set)
    out["not_ints"] = list(everything_set - ints_set)
    out["not_trues"] = list(everything_set - trues_set)
    out["not_falses"] = list(everything_set - falses_set)
    out["not_bools"] = list(everything_set - trues_set - falses_set)
    out["not_types"] = list(everything_set - types_set)
    out["not_strings"] = list(everything_set - strings_set)
    out["not_nones"] = list(everything_set - nones_set)
    out["not_nans"] = list(everything_set - nans_set)
    out["not_coercable_to_float"] = list(everything_set - coercable_to_float_set)
    out["not_coercable_to_int"] = list(everything_set - coercable_to_int_set)

    out["not_floats_or_coercable"] = list(everything_set - floats_set - coercable_to_float_set)
    out["not_ints_or_coercable"] = list(everything_set - ints_set - coercable_to_int_set)

    # Have to treat iterables specially because they are not hashable. They are 
    # not already in everything_set
    out["everything"] = list(everything_set)
    out["everything"].extend(iterables)

    out["not_iterables"] = list(everything_set - strings_set)
    for k in out.keys():
        if k.startswith("not_"):
            out[k].extend(iterables)

    # Specific types ([dict,float,int,bool,str,np.ndarray,set,tuple,pd.DataFrame])
    renamer = {"numpy.ndarray":"np.ndarray",
               "pandas.core.frame.DataFrame":"pd.DataFrame"}
    for t in types:
        if t is TestClass:
            continue

        type_name = f"{t}".split()[1][1:-2]
        if type_name in renamer:
            type_name = renamer[type_name]
        
        out[type_name] = [v for v in out["everything"] if issubclass(type(v),t)]
        out[f"not_{type_name}"] = [v for v in out["everything"] if not issubclass(type(v),t)]
    

    # Make lists with float or float-like iterables vs. not
    allowed = []
    not_allowed = []
    for t in out["everything"]:

        if issubclass(type(t),type):
            not_allowed.append(t)
            continue

        if hasattr(t,"__iter__"):

            if issubclass(type(t),str):
                try:
                    float(t)
                except ValueError:
                    not_allowed.append(t)
                    continue
            
            if issubclass(type(t),dict):
                not_allowed.append(t)
                continue

        else:  

            try:
                if np.isinf(t):
                    not_allowed.append(t)
                    continue
            except (TypeError,ValueError):
                pass

            try:
                if np.isnan(t):
                    not_allowed.append(t)
                    continue
            except (TypeError,ValueError):
                pass

            try:
                if issubclass(type(t),pd.NA):
                    not_allowed.append(t)
                    continue
            except (TypeError,ValueError):
                continue

        allowed.append(t)

    out["float_value_or_iter"] = allowed
    out["not_float_value_or_iter"] = not_allowed
    
    
    return out