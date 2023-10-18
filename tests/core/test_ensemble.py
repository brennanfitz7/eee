import pytest

from eee.core.ensemble import Ensemble

import numpy as np
import pandas as pd

def test_Ensemble(variable_types):

    ens = Ensemble()
    assert issubclass(type(ens),Ensemble)
    
    ens = Ensemble(gas_constant=1)
    assert ens._gas_constant == 1

    bad = variable_types["not_floats_or_coercable"][:]
    for v in variable_types["floats_or_coercable"]:
        if float(v) <= 0:
            bad.append(v)
            continue
        
        Ensemble(gas_constant=v)
    
    for v in bad:
        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            Ensemble(gas_constant=v)
        

def test_Ensemble_add_species(variable_types):
    
    ens = Ensemble()
    ens.add_species(name="test",
                    observable=False,
                    folded=True,
                    dG0=-1)
    
    assert len(ens._species_dict["test"]) == 3
    assert ens._species_dict["test"]["observable"] == False
    assert ens._species_dict["test"]["folded"] == True
    assert ens._species_dict["test"]["dG0"] == -1
    
    assert np.array_equal(ens._ligand_list,[])

    # Check for something already in ensemble
    with pytest.raises(ValueError):
        ens.add_species(name="test",
                        observable=False,
                        folded=True,
                        dG0=-1)
        
    ens.add_species(name="another",
                    observable=True,
                    folded=False,
                    X=1)
    
    assert ens._species_dict["another"]["observable"] == True
    assert ens._species_dict["another"]["folded"] == False
    assert ens._species_dict["another"]["dG0"] == 0
    assert ens._species_dict["another"]["X"] == 1
    assert np.array_equal(ens._ligand_list,["X"])

    # observable argument type checking
    print("--- observable ---")
    for v in variable_types["bools"]:
        print(v,type(v),flush=True)
        ens = Ensemble()
        ens.add_species(name="test",observable=v)
        assert ens._species_dict["test"]["observable"] == bool(v)

    for v in variable_types["not_bools"]:
        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            ens.add_species(name="test",observable=v)

    # folded argument type checking
    print("--- folded ---")
    for v in variable_types["bools"]:
        print(v,type(v),flush=True)
        ens = Ensemble()
        ens.add_species(name="test",folded=v)
        assert ens._species_dict["test"]["folded"] == bool(v)

    for v in variable_types["not_bools"]:
        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            ens.add_species(name="test",folded=v)

    # dG0 argument type checking
    print("--- dG0 ---")
    for v in variable_types["floats_or_coercable"]:
        print(v,type(v),flush=True)
        ens = Ensemble()
        ens.add_species(name="test",dG0=v)
        assert ens._species_dict["test"]["dG0"] == float(v)

    for v in variable_types["not_floats_or_coercable"]:
        print(v,type(v),flush=True)
        ens = Ensemble()
        with pytest.raises(ValueError):
            ens.add_species(name="test",dG0=v)

    # Pass in a ligand
    print("--- ligands ---")
    bad_ones = []
    for v in variable_types["floats_or_coercable"]:

        if float(v) < 0:
            bad_ones.append(v)
            continue

        print(v,type(v),flush=True)

        ens = Ensemble()
        ens.add_species(name="test",X=v)
        assert ens._species_dict["test"]["X"] == float(v)

    bad_ones.extend(variable_types["not_floats_or_coercable"])
    for v in bad_ones:
        print(v,type(v),flush=True)

        ens = Ensemble()
        with pytest.raises(ValueError):
            ens.add_species(name="test",X=v)


def test_Ensemble__build_z_matrix():

    # Single species, not observable, dG = 0, not coupled to ligand
    ens = Ensemble()
    ens.add_species(name="test1",
                    observable=False,
                    folded=True,
                    dG0=0)

    ens._build_z_matrix(ligand_dict={})
    assert np.array_equal(ens._z_matrix.shape,(1,1))
    assert np.array_equal(ens._z_matrix,[[0.0]])
    assert np.array_equal(ens._obs_mask,[False])
    assert np.array_equal(ens._not_obs_mask,[True])
    assert np.array_equal(ens._folded_mask,[True])
    assert np.array_equal(ens._unfolded_mask,[False])

    # Single species, not observable, dG = 0, coupled to ligand
    ens = Ensemble()
    ens.add_species(name="test1",
                    observable=False,
                    folded=True,
                    dG0=0,
                    X=1)

    ens._build_z_matrix(ligand_dict={})
    assert np.array_equal(ens._z_matrix.shape,(1,1))
    assert np.array_equal(ens._z_matrix,[[0.0]])
    assert np.array_equal(ens._obs_mask,[False])
    assert np.array_equal(ens._not_obs_mask,[True])
    assert np.array_equal(ens._folded_mask,[True])
    assert np.array_equal(ens._unfolded_mask,[False])

    # Single species, not observable, dG = 0, coupled to ligand. Now add ligand
    ens = Ensemble()
    ens.add_species(name="test1",
                    observable=False,
                    folded=True,
                    dG0=0,
                    X=1)

    ens._build_z_matrix(ligand_dict={"X":np.array([1.0])})
    assert np.array_equal(ens._z_matrix.shape,(1,1))
    assert np.array_equal(ens._z_matrix,[[-1]])
    assert np.array_equal(ens._obs_mask,[False])
    assert np.array_equal(ens._not_obs_mask,[True])
    assert np.array_equal(ens._folded_mask,[True])
    assert np.array_equal(ens._unfolded_mask,[False])

    # Two species, one observable, dG = 0, One coupled to ligand, with single ligand
    ens = Ensemble()
    ens.add_species(name="test1",
                    observable=False,
                    folded=True,
                    dG0=0,
                    X=1)
    ens.add_species(name="test2",
                    observable=True,
                    folded=False,
                    dG0=0)

    ens._build_z_matrix(ligand_dict={"X":np.array([1.0])})
    assert np.array_equal(ens._z_matrix.shape,(2,1))
    assert np.array_equal(ens._z_matrix,[[-1],[0]])
    assert np.array_equal(ens._obs_mask,[False,True])
    assert np.array_equal(ens._not_obs_mask,[True,False])
    assert np.array_equal(ens._folded_mask,[True,False])
    assert np.array_equal(ens._unfolded_mask,[False,True])

    # Two species, one observable, dG = 0, One coupled to ligand, with three ligand
    ens = Ensemble()
    ens.add_species(name="test1",
                    observable=False,
                    folded=True,
                    dG0=0,
                    X=1)
    ens.add_species(name="test2",
                    observable=True,
                    folded=False,
                    dG0=0)

    ens._build_z_matrix(ligand_dict={"X":np.array([0,0.5,1.0])})
    assert np.array_equal(ens._z_matrix.shape,(2,3))
    assert np.array_equal(ens._z_matrix,[[0,-0.5,-1],[0,0,0]])
    assert np.array_equal(ens._obs_mask,[False,True])
    assert np.array_equal(ens._not_obs_mask,[True,False])
    assert np.array_equal(ens._folded_mask,[True,False])
    assert np.array_equal(ens._unfolded_mask,[False,True])
    
    # Two species, one observable, dG = 0, 1, One coupled to ligand, with three ligand
    ens = Ensemble()
    ens.add_species(name="test1",
                    observable=False,
                    folded=True,
                    dG0=0,
                    X=1)
    ens.add_species(name="test2",
                    observable=True,
                    folded=False,
                    dG0=1)

    ens._build_z_matrix(ligand_dict={"X":np.array([0,0.5,1.0])})
    assert np.array_equal(ens._z_matrix.shape,(2,3))
    assert np.array_equal(ens._z_matrix,[[0,-0.5,-1],[1,1,1]])
    assert np.array_equal(ens._obs_mask,[False,True])
    assert np.array_equal(ens._not_obs_mask,[True,False])
    assert np.array_equal(ens._folded_mask,[True,False])
    assert np.array_equal(ens._unfolded_mask,[False,True])

    # Two species, one observable, dG = 0, 1, Both coupled to different ligand, 
    # with three ligand
    ens = Ensemble()
    ens.add_species(name="test1",
                    observable=False,
                    folded=True,
                    dG0=0,
                    X=1)
    ens.add_species(name="test2",
                    observable=True,
                    folded=False,
                    dG0=1,
                    Y=2)

    ens._build_z_matrix(ligand_dict={"X":np.array([0,0.5,1.0]),
                                 "Y":np.array([1,0.5,0.0])})
    assert np.array_equal(ens._z_matrix.shape,(2,3))
    assert np.array_equal(ens._z_matrix,[[0,-0.5,-1],[-2 + 1,-1 + 1,0 + 1]])
    assert np.array_equal(ens._obs_mask,[False,True])
    assert np.array_equal(ens._not_obs_mask,[True,False])
    assert np.array_equal(ens._folded_mask,[True,False])
    assert np.array_equal(ens._unfolded_mask,[False,True])


    # Three species, dG = 0, 1, 3 Two coupled to different ligand,  with three ligand
    ens = Ensemble()
    ens.add_species(name="test1",
                    observable=False,
                    folded=True,
                    dG0=0,
                    X=1)
    ens.add_species(name="test2",
                    observable=True,
                    folded=False,
                    dG0=1,
                    Y=2)
    ens.add_species(name="test3",
                    observable=True,
                    folded=True,
                    dG0=3)

    ens._build_z_matrix(ligand_dict={"X":np.array([0,0.5,1.0]),
                                 "Y":np.array([1,0.5,0.0])})
    assert np.array_equal(ens._z_matrix.shape,(3,3))
    assert np.array_equal(ens._z_matrix,[[0,-0.5,-1],
                                         [-2 + 1,-1 + 1,0 + 1],
                                         [3,3,3]])
    assert np.array_equal(ens._obs_mask,[False,True,True])
    assert np.array_equal(ens._not_obs_mask,[True,False,False])
    assert np.array_equal(ens._folded_mask,[True,False,True])
    assert np.array_equal(ens._unfolded_mask,[False,True,False])

def test_Ensemble__get_weights():

    # single species, R = 1, temperature = 1, no mutations
    ens = Ensemble(gas_constant=1)
    ens.add_species(name="test1",
                    observable=False,
                    folded=True,
                    dG0=0)
    
    mut_energy = np.array([0.0])
    temperature = np.ones(1,dtype=float)    

    ens._build_z_matrix(ligand_dict={})
    weights = ens._get_weights(mut_energy=mut_energy,temperature=temperature)
    assert np.array_equal(weights.shape,[1,1])
    assert np.array_equal(weights,[np.exp([ens._max_allowed])])

    # single species, R = 1, temperature = 500, no mutations. This tests the shift to 
    # max_allowed. 
    ens = Ensemble(gas_constant=1)
    ens.add_species(name="test1",
                    observable=False,
                    folded=True,
                    dG0=0)
    
    mut_energy = np.array([0.0])
    temperature = 500*np.ones(1,dtype=float)    
    ens._build_z_matrix(ligand_dict={})
    weights = ens._get_weights(mut_energy=mut_energy,temperature=temperature)
    assert np.array_equal(weights.shape,[1,1])
    assert np.array_equal(weights,[np.exp([ens._max_allowed])])
    
    # Two species. 
    ens = Ensemble(gas_constant=1)
    ens.add_species(name="test1",
                    observable=False,
                    folded=True,
                    dG0=0)
    ens.add_species(name="test2",
                    observable=False,
                    folded=False,
                    dG0=0)
    
    mut_energy = np.array([0.0,0.0])
    temperature = np.ones(1,dtype=float)    
    ens._build_z_matrix(ligand_dict={})
    weights = ens._get_weights(mut_energy=mut_energy,temperature=temperature)
    assert np.array_equal(weights.shape,[2,1])
    assert np.array_equal(weights,[[np.exp(ens._max_allowed)],
                                   [np.exp(ens._max_allowed)]])


    # Two species. Mutate one. 
    ens = Ensemble(gas_constant=1)
    ens.add_species(name="test1",
                    observable=False,
                    folded=True,
                    dG0=0)
    ens.add_species(name="test2",
                    observable=False,
                    folded=False,
                    dG0=0)
    
    mut_energy = np.array([0.0,1.0])
    temperature = np.ones(1,dtype=float)    
    ens._build_z_matrix(ligand_dict={})
    weights = ens._get_weights(mut_energy=mut_energy,temperature=temperature)
    assert np.array_equal(weights.shape,[2,1])
    norm_weights = weights/np.sum(weights,axis=0)
    Z = np.exp(0) + np.exp(-1)
    assert np.isclose(norm_weights[0,0],np.exp(0)/Z)
    assert np.isclose(norm_weights[1,0],np.exp(-1)/Z)

    # Make sure temperature works as expected
    temperature = 50*np.ones(1,dtype=float)                                        
    weights = ens._get_weights(mut_energy=mut_energy,temperature=temperature)
    assert np.array_equal(weights.shape,[2,1])
    norm_weights = weights/np.sum(weights,axis=0)
    Z = np.exp(0/50) + np.exp(-1/50)
    assert np.isclose(norm_weights[0,0],np.exp(0/50)/Z)
    assert np.isclose(norm_weights[1,0],np.exp(-1/50)/Z)
                                        

    # Two species. Mutate one. Alter R to test gas constant
    ens = Ensemble(gas_constant=50)
    ens.add_species(name="test1",
                    observable=False,
                    folded=True,
                    dG0=0)
    ens.add_species(name="test2",
                    observable=False,
                    folded=False,
                    dG0=0)
    
    mut_energy = np.array([0.0,1.0])
    temperature = np.ones(1,dtype=float)*50    
    ens._build_z_matrix(ligand_dict={})
                                  
    weights = ens._get_weights(mut_energy=mut_energy,temperature=temperature)
    assert np.array_equal(weights.shape,[2,1])
    norm_weights = weights/np.sum(weights,axis=0)
    Z = np.exp(0/2500) + np.exp(-1/2500)
    assert np.isclose(norm_weights[0,0],np.exp(0/2500)/Z)
    assert np.isclose(norm_weights[1,0],np.exp(-1/2500)/Z)


    # Two species. dG0
    ens = Ensemble(gas_constant=1)
    ens.add_species(name="test1",
                    observable=False,
                    folded=True,
                    dG0=-1)
    ens.add_species(name="test2",
                    observable=False,
                    folded=False,
                    dG0=0)
    
    mut_energy = np.array([0.0,1.0])
    temperature = np.ones(1,dtype=float)    
    ens._build_z_matrix(ligand_dict={})
                        
    weights = ens._get_weights(mut_energy=mut_energy,temperature=temperature)
    assert np.array_equal(weights.shape,[2,1])
    norm_weights = weights/np.sum(weights,axis=0)
    Z = np.exp(1) + np.exp(-1)
    assert np.isclose(norm_weights[0,0],np.exp(1)/Z)
    assert np.isclose(norm_weights[1,0],np.exp(-1)/Z)

    # Two species. dG0. Add ligand_dict perturbation
    ens = Ensemble(gas_constant=1)
    ens.add_species(name="test1",
                    observable=False,
                    folded=True,
                    dG0=-1,
                    X=1)
    ens.add_species(name="test2",
                    observable=False,
                    folded=False,
                    dG0=0)
    
    mut_energy = np.array([0.0,1.0])
    temperature = np.ones(1,dtype=float)    
    ens._build_z_matrix(ligand_dict={"X":np.array([0,0.5,1])})
                        
    weights = ens._get_weights(mut_energy=mut_energy,temperature=temperature)
    assert np.array_equal(weights.shape,[2,3])
    norm_weights = weights/np.sum(weights,axis=0)
    Z = np.exp(1) + np.exp(-1)
    assert np.isclose(norm_weights[0,0],np.exp(1)/Z)
    assert np.isclose(norm_weights[1,0],np.exp(-1)/Z)

    Z = np.exp(1 + 0.5) + np.exp(-1)
    assert np.isclose(norm_weights[0,1],np.exp(1+0.5)/Z)
    assert np.isclose(norm_weights[1,1],np.exp(-1)/Z)

    Z = np.exp(1 + 1.0) + np.exp(-1)
    assert np.isclose(norm_weights[0,2],np.exp(1+1.0)/Z)
    assert np.isclose(norm_weights[1,2],np.exp(-1)/Z)

def test_Ensemble_get_species_dG(variable_types):

    ens = Ensemble()
    ens.add_species(name="test",
                    observable=False,
                    folded=True,
                    dG0=-1)

    with pytest.raises(ValueError):
        ens.get_species_dG("not_a_species")
    
    assert ens.get_species_dG("test") == -1
    assert ens.get_species_dG("test",mut_energy=10) == 9
    assert ens.get_species_dG("test",mut_energy=10,ligand_dict={"X":10}) == 9

    ens.add_species(name="another",
                    observable=True,
                    folded=False,
                    X=1)
    
    assert ens.get_species_dG("another") == 0
    assert ens.get_species_dG("another",mut_energy=10) == 10
    assert ens.get_species_dG("another",mut_energy=10,ligand_dict={"X":10}) == 0

    # Pass in an array of ligand X
    value = ens.get_species_dG("another",mut_energy=10,ligand_dict={"X":np.arange(10)})
    assert np.array_equal(value,-np.arange(10) + 10)

    # Stoichiometry of 2
    ens = Ensemble()
    ens.add_species(name="stoich2",
                    observable=False,
                    folded=True,
                    dG0=-5,
                    X=2)
    value = ens.get_species_dG("stoich2",mut_energy=10,ligand_dict={"X":np.arange(10)})
    assert np.array_equal(value,-np.arange(10)*2 + 5)

    # Two chemical potentials same species
    ens = Ensemble()
    ens.add_species(name="two_ligand",
                    observable=False,
                    folded=True,
                    dG0=-5,
                    X=2,
                    Y=1)
    value = ens.get_species_dG("two_ligand",
                               mut_energy=0,
                               ligand_dict={"X":np.arange(10),
                                        "Y":np.arange(10)})
    assert np.array_equal(value,-np.arange(10)*2+-np.arange(10)-5)

    # mut_energy argument type checking
    print("--- mut_energy ---")
    for v in variable_types["floats_or_coercable"]:
        print(v,type(v),flush=True)
        ens = Ensemble()
        ens.add_species(name="test")
        assert ens.get_species_dG(name="test",mut_energy=v) == float(v)

    for v in variable_types["not_floats_or_coercable"]:
        print(v,type(v),flush=True)
        ens = Ensemble()
        ens.add_species(name="test")
        with pytest.raises(ValueError):
            ens.get_species_dG(name="test",mut_energy=v)

    # ligand_dict argument type checking
    print("--- ligand_dict ---")
    for v in variable_types["dict"]:
        print(v,type(v),flush=True)
        ens = Ensemble()
        ens.add_species(name="test1")
        ens.add_species(name="test2",observable=True)
        ens.get_obs(ligand_dict=v)

    for v in variable_types["not_dict"]:
        print(v,type(v),flush=True)

        # none okay
        if v is None:
            continue

        ens = Ensemble()
        ens.add_species(name="test1")
        ens.add_species(name="test2",observable=True)
        with pytest.raises(ValueError):
            ens.get_obs(ligand_dict=v)

    for v in variable_types["float_value_or_iter"]:
        print(v,type(v),flush=True)
        
        if hasattr(v,"__iter__") and len(v) == 0:
            continue

        if issubclass(type(v),pd.DataFrame):
            continue

        ligand_dict = {"X":v}
        ens = Ensemble()
        ens.add_species(name="test1")
        ens.add_species(name="test2",observable=True)
        ens.get_obs(ligand_dict=ligand_dict)

    not_allowed = variable_types["not_float_value_or_iter"][:]
    not_allowed.append([])
    not_allowed.append(pd.DataFrame({"X":[1,2,3]}))

    for v in not_allowed:
        print(v,type(v),flush=True)
    
        ligand_dict = {"X":v}
        ens = Ensemble()
        ens.add_species(name="test1")
        ens.add_species(name="test2",observable=True)
        with pytest.raises(ValueError):
            ens.get_obs(ligand_dict=ligand_dict)


def test_Ensemble_get_obs(variable_types):

    # ------------------------------------------------------------------------
    # Not enough species
    ens = Ensemble()
    ens.add_species(name="test",
                    folded=True,
                    observable=False)
    with pytest.raises(ValueError):
        ens.get_obs()

    # ------------------------------------------------------------------------
    # No observable species
    ens = Ensemble()
    ens.add_species(name="test1",
                    folded=True,
                    observable=False)
    ens.add_species(name="test2",
                    folded=False,
                    observable=False)
    with pytest.raises(ValueError):
        ens.get_obs()

    # ------------------------------------------------------------------------
    # Only observable species
    ens = Ensemble()
    ens.add_species(name="test1",
                    observable=True,
                    folded=True)
    ens.add_species(name="test2",
                    observable=True,
                    folded=False)
    with pytest.raises(ValueError):
        ens.get_obs()

    # ------------------------------------------------------------------------
    # One observable, one not
    ens = Ensemble()
    ens.add_species(name="test1",
                    observable=True,
                    folded=True)
    ens.add_species(name="test2",
                    observable=False,
                    folded=False)
    
    df = ens.get_obs()
    assert df.loc[0,"dG_obs"] == 0
    assert df.loc[0,"fx_obs"] == 0.5
    assert df.loc[0,"test1"] == 0.5
    assert df.loc[0,"test2"] == 0.5
    assert df.loc[0,"fx_folded"] == 0.5

    # ------------------------------------------------------------------------
    # One observable, two not. (Use R = 1 and temperature = 1 to simplify math)
    ens = Ensemble(gas_constant=1)
    ens.add_species(name="test1",
                    observable=True,
                    folded=False)
    ens.add_species(name="test2",
                    observable=False,
                    folded=True)
    ens.add_species(name="test3",   
                    observable=False,
                    folded=True)
    
    df = ens.get_obs(temperature=1)
    assert df.loc[0,"dG_obs"] == -np.log(1/2)
    assert df.loc[0,"fx_obs"] == 1/3
    assert df.loc[0,"test1"] == 1/3
    assert df.loc[0,"test2"] == 1/3
    assert df.loc[0,"fx_folded"] == 2/3

    # ------------------------------------------------------------------------
    # One observable, two not. (Use R = 1 and temperature = 1 to simplify math). ligand_dict
    # interesting.
    ens = Ensemble(gas_constant=1)
    ens.add_species(name="test1",
                    observable=True,
                    folded=False,
                    X=1)
    ens.add_species(name="test2",
                    observable=False,
                    folded=True)
    ens.add_species(name="test3",   
                    observable=False,
                    folded=True)
    df = ens.get_obs(temperature=1,
                     ligand_dict={"X":[0,1]})
    
    test1 = np.exp(0)
    test2 = np.exp(0)
    test3 = np.exp(0)
    test_all = test1 + test2 + test3
    numerator = test1
    denominator = test2 + test3

    dG = -np.log(numerator/denominator)
    fx = numerator/(numerator + denominator)

    assert np.isclose(df.loc[0,"dG_obs"],dG)
    assert np.isclose(df.loc[0,"fx_obs"],fx)
    assert np.isclose(df.loc[0,"test1"],test1/test_all)
    assert np.isclose(df.loc[0,"test2"],test2/test_all)
    assert np.isclose(df.loc[0,"test3"],test3/test_all)
    assert np.isclose(df.loc[0,"fx_folded"],1-fx)

    test1 = np.exp(1)
    test2 = np.exp(0)
    test3 = np.exp(0)
    test_all = test1 + test2 + test3
    numerator = test1
    denominator = test2 + test3

    dG = -np.log(numerator/denominator)
    fx = numerator/(numerator + denominator)

    assert np.isclose(df.loc[1,"dG_obs"],dG)
    assert np.isclose(df.loc[1,"fx_obs"],fx)
    assert np.isclose(df.loc[1,"test1"],test1/test_all)
    assert np.isclose(df.loc[1,"test2"],test2/test_all)
    assert np.isclose(df.loc[1,"test3"],test3/test_all)
    assert np.isclose(df.loc[1,"fx_folded"],1-fx)

    # ------------------------------------------------------------------------
    # One observable, two not. (Use R = 1 and temperature = 1 to simplify math). ligand_dict
    # interesting. mut_energy interesting
    ens = Ensemble(gas_constant=1)
    ens.add_species(name="test1",
                    observable=True,
                    folded=False,
                    X=1)
    ens.add_species(name="test2",
                    observable=False,
                    folded=True)
    ens.add_species(name="test3",   
                    observable=False,
                    folded=True)
    df = ens.get_obs(temperature=1,
                     ligand_dict={"X":[0,1]},
                     mut_energy={"test2":-1})
    
    test1 = np.exp(0)
    test2 = np.exp(1)
    test3 = np.exp(0)
    test_all = test1 + test2 + test3
    numerator = test1
    denominator = test2 + test3

    dG = -np.log(numerator/denominator)
    fx = numerator/(numerator + denominator)

    assert np.isclose(df.loc[0,"dG_obs"],dG)
    assert np.isclose(df.loc[0,"fx_obs"],fx)
    assert np.isclose(df.loc[0,"test1"],test1/test_all)
    assert np.isclose(df.loc[0,"test2"],test2/test_all)
    assert np.isclose(df.loc[0,"test3"],test3/test_all)
    assert np.isclose(df.loc[0,"fx_folded"],1-fx)

    test1 = np.exp(1)
    test2 = np.exp(1)
    test3 = np.exp(0)
    test_all = test1 + test2 + test3
    numerator = test1
    denominator = test2 + test3

    dG = -np.log(numerator/denominator)
    fx = numerator/(numerator + denominator)

    assert np.isclose(df.loc[1,"dG_obs"],dG)
    assert np.isclose(df.loc[1,"fx_obs"],fx)
    assert np.isclose(df.loc[1,"test1"],test1/test_all)
    assert np.isclose(df.loc[1,"test2"],test2/test_all)
    assert np.isclose(df.loc[1,"test3"],test3/test_all)
    assert np.isclose(df.loc[1,"fx_folded"],1-fx)
    
    # ------------------------------------------------------------------------
    # One observable, two not. (Use R = 1 and temperature = 1 to simplify math). ligand_dict
    # interesting, diff for different species. mut_energy interesting
    ens = Ensemble(gas_constant=1)
    ens.add_species(name="test1",
                    observable=True,
                    folded=True,
                    X=1)
    ens.add_species(name="test2",
                    observable=False,
                    folded=False,
                    Y=1)
    ens.add_species(name="test3",   
                    observable=False,
                    folded=False)
    df = ens.get_obs(temperature=1,
                     ligand_dict={"X":[0,1],"Y":[0,1]},
                     mut_energy={"test2":-1})
    
    test1 = np.exp(0)
    test2 = np.exp(1)
    test3 = np.exp(0)
    test_all = test1 + test2 + test3
    numerator = test1
    denominator = test2 + test3

    dG = -np.log(numerator/denominator)
    fx = numerator/(numerator + denominator)

    assert np.isclose(df.loc[0,"dG_obs"],dG)
    assert np.isclose(df.loc[0,"fx_obs"],fx)
    assert np.isclose(df.loc[0,"test1"],test1/test_all)
    assert np.isclose(df.loc[0,"test2"],test2/test_all)
    assert np.isclose(df.loc[0,"test3"],test3/test_all)
    assert np.isclose(df.loc[0,"fx_folded"],fx)

    test1 = np.exp(1)
    test2 = np.exp(2)
    test3 = np.exp(0)
    test_all = test1 + test2 + test3
    numerator = test1
    denominator = test2 + test3

    dG = -np.log(numerator/denominator)
    fx = numerator/(numerator + denominator)

    assert np.isclose(df.loc[1,"dG_obs"],dG)
    assert np.isclose(df.loc[1,"fx_obs"],fx)
    assert np.isclose(df.loc[1,"test1"],test1/test_all)
    assert np.isclose(df.loc[1,"test2"],test2/test_all)
    assert np.isclose(df.loc[1,"test3"],test3/test_all)
    assert np.isclose(df.loc[1,"fx_folded"],fx)


    # ------------------------------------------------------------------------
    # Overflow protection

    # this dG will overflow if put in raw but should be fine if we shift. 
    max_allowed = np.log(np.finfo("d").max)

    # Make sure this should overflow if used naively
    with pytest.warns(RuntimeWarning):
        x = np.exp(max_allowed + 1)
    assert np.isinf(x)

    ens = Ensemble(gas_constant=1)
    ens.add_species(name="test1",
                    observable=True,
                    folded=True,
                    dG0=max_allowed+2)
    ens.add_species(name="test2",
                    observable=False,
                    folded=False,
                    dG0=max_allowed+1)
    
    df = ens.get_obs(temperature=1)
    assert np.isclose(df.loc[0,"dG_obs"],-np.log(np.exp(-1)/np.exp(0)))
    assert np.isclose(df.loc[0,"fx_obs"],np.exp(-1)/(np.exp(0) + np.exp(-1)))
    assert np.isclose(df.loc[0,"test1"],np.exp(-1)/(np.exp(0) + np.exp(-1)))
    assert np.isclose( df.loc[0,"test2"],np.exp(0)/(np.exp(0) + np.exp(-1)))
    assert np.isclose(df.loc[0,"fx_folded"],np.exp(-1)/(np.exp(0) + np.exp(-1)))


    # ligand_dict argument type checking
    print("--- ligand_dict ---")
    for v in variable_types["dict"]:
        print(v,type(v),flush=True)
        ens = Ensemble()
        ens.add_species(name="test1")
        ens.add_species(name="test2",observable=True)
        ens.get_obs(ligand_dict=v)

    for v in variable_types["not_dict"]:
        print(v,type(v),flush=True)

        # none okay
        if v is None:
            continue

        ens = Ensemble()
        ens.add_species(name="test1")
        ens.add_species(name="test2",observable=True)
        with pytest.raises(ValueError):
            ens.get_obs(ligand_dict=v)

    for v in variable_types["float_value_or_iter"]:
        print(v,type(v),flush=True)
        
        if hasattr(v,"__iter__") and len(v) == 0:
            continue

        if issubclass(type(v),pd.DataFrame):
            continue

        ligand_dict = {"X":v}
        ens = Ensemble()
        ens.add_species(name="test1")
        ens.add_species(name="test2",observable=True)
        ens.get_obs(ligand_dict=ligand_dict)

    not_allowed = variable_types["not_float_value_or_iter"][:]
    not_allowed.append([])
    not_allowed.append(pd.DataFrame({"X":[1,2,3]}))

    for v in not_allowed:
        print(v,type(v),flush=True)
    
        ligand_dict = {"X":v}
        ens = Ensemble()
        ens.add_species(name="test1")
        ens.add_species(name="test2",observable=True)
        with pytest.raises(ValueError):
            ens.get_obs(ligand_dict=ligand_dict)

    print("--- mut_energy ---")
    for v in [{},{"test1":1},{"test2":1},{"test1":1,"test2":1}]: 
        print(v,type(v),flush=True)

        ens = Ensemble()
        ens.add_species(name="test1")
        ens.add_species(name="test2",observable=True)
        ens.get_obs(mut_energy=v)

    not_allowed = variable_types["not_dict"]
    for v in not_allowed:
        print(v,type(v),flush=True)
        if v is None:
            continue

        ens = Ensemble()
        ens.add_species(name="test1")
        ens.add_species(name="test2",observable=True)
        with pytest.raises(ValueError):
            ens.get_obs(mut_energy=v)

    for v in variable_types["floats_or_coercable"]:
        print(v,type(v),flush=True)
        mut_energy = {"test1":v}

        ens = Ensemble()
        ens.add_species(name="test1")
        ens.add_species(name="test2",observable=True)
        ens.get_obs(mut_energy=mut_energy)
    

    for v in variable_types["not_floats_or_coercable"]:
        print(v,type(v),flush=True)
        mut_energy = {"test1":v}

        ens = Ensemble()
        ens.add_species(name="test1")
        ens.add_species(name="test2",observable=True)
        with pytest.raises(ValueError):
            ens.get_obs(mut_energy=mut_energy)


    print("--- temperature ---")
    for v in [1,"1",1.0]:
        print(v,type(v),flush=True)

        ens = Ensemble()
        ens.add_species(name="test1")
        ens.add_species(name="test2",observable=True)
        ens.get_obs(temperature=v)

    not_allowed = variable_types["not_floats_or_coercable"][:]
    not_allowed.append(0.0)
    not_allowed.append(-1.0)
    for v in not_allowed:
        print(v,type(v),flush=True)

        ens = Ensemble()
        ens.add_species(name="test1")
        ens.add_species(name="test2",observable=True)
        with pytest.raises(ValueError):
            ens.get_obs(temperature=v)


def test_Ensemble_read_ligand_dict(variable_types):

    # Two species. dG0. Add ligand_dict perturbation
    ens = Ensemble()
    ens.add_species(name="test1",
                    observable=False,
                    folded=True,
                    dG0=0,
                    X=1)
    ens.add_species(name="test2",
                    observable=True,
                    folded=False,
                    dG0=1,
                    Y=2)
    
    for v in variable_types["not_dict"]:
        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            ens.read_ligand_dict(v)

    assert not hasattr(ens,"_z_matrix")

    # Just check one load -- wraps _create_z_matrix which we already test 
    # extensively. 
    ens.read_ligand_dict(ligand_dict={"X":np.array([0,0.5,1.0]),
                              "Y":np.array([1,0.5,0.0])})
    assert np.array_equal(ens._z_matrix.shape,(2,3))
    assert np.array_equal(ens._z_matrix,[[0,-0.5,-1],[-2 + 1,-1 + 1,0 + 1]])
    assert np.array_equal(ens._obs_mask,[False,True])
    assert np.array_equal(ens._not_obs_mask,[True,False])
    assert np.array_equal(ens._folded_mask,[True,False])
    assert np.array_equal(ens._unfolded_mask,[False,True])


def test_Ensemble_mut_dict_to_array():

    # Two species. 
    ens = Ensemble(gas_constant=1)
    ens.add_species(name="test1",
                    observable=False,
                    dG0=0)
    ens.add_species(name="test2",
                    observable=False,
                    dG0=0)
    
    out_array = ens.mut_dict_to_array({"test1":1.0,"test2":2.0})
    assert np.array_equal(out_array,[1,2])

    out_array = ens.mut_dict_to_array({"test2":1.0,"test1":2.0})
    assert np.array_equal(out_array,[2,1])

    # Hack that should invert outputs. Never really happen in real life.
    ens._species_list = ["test2","test1"]
    out_array = ens.mut_dict_to_array({"test2":1.0,"test1":2.0})
    assert np.array_equal(out_array,[1,2])

def test_Ensemble_get_fx_obs_fast():

    ens = Ensemble(gas_constant=1)
    ens.add_species(name="test1",
                    observable=False,
                    folded=True,
                    dG0=1,
                    X=1)
    ens.add_species(name="test2",
                    observable=True,
                    folded=False,
                    dG0=0)
    ens.add_species(name="test3",
                    observable=True,
                    folded=True,
                    dG0=2)

    ens.read_ligand_dict(ligand_dict={"X":np.array([0,1.0])})
    temperature = np.ones(1,dtype=float)    

    value, fx_folded = ens.get_fx_obs_fast(mut_energy_array=np.array([0,0,0]),
                                           temperature=temperature)
    t1 = np.exp(-1)
    t2 = np.exp(-0)
    t3 = np.exp(-2)
    Z = t1 + t2 + t3
    predicted = [(t2 + t3)/Z]
    t1 = np.exp(0)
    t2 = np.exp(-0)
    t3 = np.exp(-2)
    Z = t1 + t2 + t3
    predicted.append((t2 + t3)/Z)

    assert np.array_equal(np.round(value,2),
                          np.round(predicted,2))

    t1 = np.exp(-1) 
    t2 = np.exp(-0)
    t3 = np.exp(-2)
    Z = t1 + t2 + t3
    predicted = [(t1 + t3)/Z]
    t1 = np.exp(0)
    t2 = np.exp(-0)
    t3 = np.exp(-2)
    Z = t1 + t2 + t3
    predicted.append((t1 + t3)/Z)

    assert np.array_equal(np.round(fx_folded,2),
                          np.round(predicted,2))

    value, fx_folded = ens.get_fx_obs_fast(mut_energy_array=np.array([0,-1,0]),
                                           temperature=temperature)
    t1 = np.exp(-1)
    t2 = np.exp(1)
    t3 = np.exp(-2)
    Z = t1 + t2 + t3
    predicted = [(t2 + t3)/Z]
    t1 = np.exp(0)
    t2 = np.exp(1)
    t3 = np.exp(-2)
    Z = t1 + t2 + t3
    predicted.append((t2 + t3)/Z)

    assert np.array_equal(np.round(value,2),
                          np.round(predicted,2))
    
    t1 = np.exp(-1)
    t2 = np.exp(1)
    t3 = np.exp(-2)
    Z = t1 + t2 + t3
    predicted = [(t1 + t3)/Z]
    t1 = np.exp(0)
    t2 = np.exp(1)
    t3 = np.exp(-2)
    Z = t1 + t2 + t3
    predicted.append((t1 + t3)/Z)

    assert np.array_equal(np.round(fx_folded,2),
                          np.round(predicted,2))


def test_Ensemble_get_dG_obs_fast():

    
    ens = Ensemble(gas_constant=1)
    ens.add_species(name="test1",
                    observable=False,
                    folded=True,
                    dG0=1,
                    X=1)
    ens.add_species(name="test2",
                    observable=True,
                    folded=False,
                    dG0=0)
    ens.add_species(name="test3",
                    observable=True,
                    folded=True,
                    dG0=2)

    ens.read_ligand_dict(ligand_dict={"X":np.array([0,1.0])})
    temperature = np.ones(1,dtype=float)    

    value, fx_folded = ens.get_dG_obs_fast(mut_energy_array=np.array([0,0,0]),temperature=temperature)
    t1 = np.exp(-1)
    t2 = np.exp(-0)
    t3 = np.exp(-2)
    predicted = [-np.log((t2 + t3)/t1)]

    t1 = np.exp(0)
    t2 = np.exp(-0)
    t3 = np.exp(-2)
    predicted.append(-np.log((t2 + t3)/t1))

    assert np.array_equal(np.round(value,2),
                          np.round(predicted,2))
    
    t1 = np.exp(-1)
    t2 = np.exp(-0)
    t3 = np.exp(-2)
    predicted = [(t1 + t3)/(t1+t2+t3)]

    t1 = np.exp(0)
    t2 = np.exp(-0)
    t3 = np.exp(-2)
    predicted.append((t1 + t3)/(t1+t2+t3))

    assert np.array_equal(np.round(fx_folded,2),
                          np.round(predicted,2))
    

    value, fx_folded = ens.get_dG_obs_fast(mut_energy_array=np.array([0,-1,0]),temperature=temperature)
    t1 = np.exp(-1)
    t2 = np.exp(1)
    t3 = np.exp(-2)
    predicted = [-np.log((t2 + t3)/t1)]

    t1 = np.exp(0)
    t2 = np.exp(1)
    t3 = np.exp(-2)
    predicted.append(-np.log((t2 + t3)/t1))

    assert np.array_equal(np.round(value,2),
                          np.round(predicted,2))

    t1 = np.exp(-1)
    t2 = np.exp(1)
    t3 = np.exp(-2)
    predicted = [(t1 + t3)/(t1+t2+t3)]

    t1 = np.exp(0)
    t2 = np.exp(1)
    t3 = np.exp(-2)
    predicted.append((t1 + t3)/(t1+t2+t3))

    assert np.array_equal(np.round(fx_folded,2),
                          np.round(predicted,2))

def test_Ensemble_to_dict():

    ens = Ensemble(gas_constant=1)
    ens.add_species(name="test1",
                    observable=False,
                    folded=True,
                    dG0=1,
                    X=1)
    
    out = ens.to_dict()
    assert out["gas_constant"] == 1
    assert len(out) == 2
    assert out["ens"]["test1"]["observable"] == False
    assert out["ens"]["test1"]["folded"] == True
    assert out["ens"]["test1"]["dG0"] == 1
    assert out["ens"]["test1"]["X"] == 1
    assert len(out["ens"]) == 1
    assert len(out["ens"]["test1"]) == 4

def test_Ensemble_get_observable_function(variable_types):

    ens = Ensemble(gas_constant=1)
    assert ens.get_observable_function("fx_obs") == ens.get_fx_obs_fast
    assert ens.get_observable_function("dG_obs") == ens.get_dG_obs_fast

    for v in variable_types["everything"]:
        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            ens.get_observable_function(v)

def test_Ensemble_species():

    ens = Ensemble(gas_constant=1)
    assert len(ens.species) == 0
    ens.add_species("test1")
    ens.add_species("test2")

    assert np.array_equal(ens.species,["test1","test2"])

def test_Ensemble_ligands():
    
    ens = Ensemble(gas_constant=1)
    assert len(ens.ligands) == 0
    ens.add_species("test1",X=1)
    ens.add_species("test2",Y=1)

    assert np.array_equal(ens.ligands,["X","Y"])

def test_Ensemble_species_df():

    ens = Ensemble(gas_constant=1)
    assert len(ens.species_df) == 0
    assert issubclass(type(ens.species_df),pd.DataFrame)

    ens.add_species("test1",X=1,dG0=5,observable=True)
    ens.add_species("test2",Y=1,folded=False,observable=False)

    assert np.array_equal(ens.species_df["name"],["test1","test2"])
    assert np.array_equal(ens.species_df["dG0"],[5,0])
    assert np.array_equal(ens.species_df["folded"],[True,False])
    assert np.array_equal(ens.species_df["observable"],[True,False])
    assert np.array_equal(ens.species_df["X"],[1,0])
    assert np.array_equal(ens.species_df["Y"],[0,1])

