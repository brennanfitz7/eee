
from eee.core.genotype.single_genotype import SingleGenotype

import numpy as np


def test_SingleGenotype(ens_test_data):

    ens = ens_test_data["ens"]
    ddg_dict = ens_test_data["ddg_dict"]

    # ---------------------------------
    # minimal input
    g = SingleGenotype(ens=ens,
                 ddg_dict=ddg_dict)
    
    assert g._ens is ens
    assert g._ddg_dict is ddg_dict
    assert len(g.sites) == 0
    assert len(g.mutations) == 0
    assert len(g.mutations_accumulated) == 0
    assert len(g.mut_energy) == 2
    assert g.mut_energy[0] == 0
    assert g.mut_energy[1] == 0

    # ---------------------------------
    # sites, mutations, mut_energy, and mutations_accumulated input
    g = SingleGenotype(ens=ens,
                 ddg_dict=ddg_dict,
                 sites=[1],
                 mutations=["M1A"],
                 mutations_accumulated=["M1A"],
                 mut_energy=np.array([1,-1]))
    
    assert g._ens is ens
    assert g._ddg_dict is ddg_dict
    assert len(g.sites) == 1
    assert g.sites[0] == 1
    assert len(g.mutations) == 1
    assert g.mutations[0] == "M1A"
    assert len(g.mutations_accumulated) == 1
    assert g.mutations_accumulated[0] == "M1A"
    assert len(g.mut_energy) == 2
    assert g.mut_energy[0] == 1
    assert g.mut_energy[1] == -1


def test_SingleGenotype_copy(ens_test_data):
    
    ens = ens_test_data["ens"]
    ddg_dict = ens_test_data["ddg_dict"]

    # ---------------------------------
    # Create a genotype object
    g = SingleGenotype(ens=ens,
                 ddg_dict=ddg_dict,
                 sites=[1],
                 mutations=["M1A"],
                 mutations_accumulated=["M1A"],
                 mut_energy=np.array([1,-1]))
    
    assert g._ens is ens
    assert g._ddg_dict is ddg_dict
    assert len(g.sites) == 1
    assert g.sites[0] == 1
    assert len(g.mutations) == 1
    assert g.mutations[0] == "M1A"
    assert len(g.mutations_accumulated) == 1
    assert g.mutations_accumulated[0] == "M1A"
    assert len(g.mut_energy) == 2
    assert g.mut_energy[0] == 1
    assert g.mut_energy[1] == -1

    # ---------------------------------
    # Now test copying of that genotype

    g2 = g.copy()

    # These should be the same object (copy by reference)
    assert g._ens is g2._ens
    assert g._ddg_dict is g2._ddg_dict

    # These should be different objects
    assert g._sites is not g2._sites
    assert g._mutations is not g2._mutations
    assert g._mutations_accumulated is not g2._mutations_accumulated
    assert g._mut_energy is not g2._mut_energy

    # But they should have copied contents correctly
    assert len(g2.sites) == 1
    assert g2.sites[0] == 1
    assert len(g2.mutations) == 1
    assert g2.mutations[0] == "M1A"
    assert len(g2.mutations_accumulated) == 1
    assert g2.mutations_accumulated[0] == "M1A"
    assert len(g2.mut_energy) == 2
    assert g2.mut_energy[0] == 1
    assert g2.mut_energy[1] == -1

def test_SingleGenotype_mutate(ens_test_data):
    
    ens = ens_test_data["ens"]
    ddg_dict = ens_test_data["ddg_dict"]

    # ---------------------------------
    # Create a genotype object
    g = SingleGenotype(ens=ens,
                 ddg_dict=ddg_dict,
                 sites=[1],
                 mutations=["M1A"],
                 mutations_accumulated=["M1A"],
                 mut_energy=np.array([1,-1]))
    
    assert len(g.sites) == 1
    assert g.sites[0] == 1
    assert len(g.mutations) == 1
    assert g.mutations[0] == "M1A"
    assert len(g.mutations_accumulated) == 1
    assert g.mutations_accumulated[0] == "M1A"
    assert len(g.mut_energy) == 2
    assert g.mut_energy[0] == 1
    assert g.mut_energy[1] == -1

    # Copy and mutate new site, new mutation
    g2 = g.copy()
    g2.mutate(2,"P2R")
    assert len(g2.sites) == 2
    assert g2.sites[0] == 1
    assert g2.sites[1] == 2
    assert len(g2.mutations) == 2
    assert g2.mutations[0] == "M1A"
    assert g2.mutations[1] == "P2R"
    assert len(g2.mutations_accumulated) == 2
    assert g2.mutations_accumulated[0] == "M1A"
    assert g2.mutations_accumulated[1] == "P2R"
    assert len(g2.mut_energy) == 2
    assert g2.mut_energy[0] == 1
    assert g2.mut_energy[1] == 0 # 1 unit higher 

    # Copy and mutate same site, new mutation
    g3 = g.copy()
    g3.mutate(1,"M1V")
    assert len(g3.sites) == 1
    assert g3.sites[0] == 1
    assert len(g3.mutations) == 1
    assert g3.mutations[0] == "M1V"
    assert len(g3.mutations_accumulated) == 2
    assert g3.mutations_accumulated[0] == "M1A"
    assert g3.mutations_accumulated[1] == "A1V"
    assert len(g3.mut_energy) == 2
    assert g3.mut_energy[0] == -1
    assert g3.mut_energy[1] == 1

    # Copy and mutate same site, same mutation --> reversion
    g4 = g.copy()
    g4.mutate(1,"M1A")
    assert len(g4.sites) == 0
    assert len(g4.mutations) == 0
    assert len(g4.mutations_accumulated) == 2
    assert g4.mutations_accumulated[0] == "M1A"
    assert g4.mutations_accumulated[1] == "A1M"
    assert len(g4.mut_energy) == 2
    assert g4.mut_energy[0] == 0
    assert g4.mut_energy[1] == 0

    # Start at wildtype and mutate
    g = SingleGenotype(ens=ens,
                 ddg_dict=ddg_dict)
    
    g.mutate(2,"P2R")
    assert len(g.sites) == 1
    assert g.sites[0] == 2
    assert len(g.mutations) == 1
    assert g.mutations[0] == "P2R"
    assert len(g.mutations_accumulated) == 1
    assert g.mutations_accumulated[0] == "P2R"
    assert len(g.mut_energy) == 2
    assert g.mut_energy[0] == 0
    assert g.mut_energy[1] == 1

def test_SingleGenotype_mut_energy(ens_test_data):
    
    ens = ens_test_data["ens"]
    ddg_dict = ens_test_data["ddg_dict"]

    # ---------------------------------
    # Create a genotype object
    g = SingleGenotype(ens=ens,
                 ddg_dict=ddg_dict,
                 sites=[1],
                 mutations=["M1A"],
                 mutations_accumulated=["M1A"],
                 mut_energy=np.array([1,-1]))
    
    assert g.mut_energy is g._mut_energy
    assert g.mut_energy[0] == 1
    assert g.mut_energy[1] == -1

def test_SingleGenotype_mutations(ens_test_data):
    
    ens = ens_test_data["ens"]
    ddg_dict = ens_test_data["ddg_dict"]

    # ---------------------------------
    # Create a genotype object
    g = SingleGenotype(ens=ens,
                 ddg_dict=ddg_dict,
                 sites=[1],
                 mutations=["M1A"],
                 mutations_accumulated=["M1A"],
                 mut_energy=np.array([1,-1]))
    
    assert g.mutations is g._mutations
    assert np.array_equal(g.mutations,["M1A"])


def test_SingleGenotype_mutations_accumulated(ens_test_data):
    
    ens = ens_test_data["ens"]
    ddg_dict = ens_test_data["ddg_dict"]

    # ---------------------------------
    # Create a genotype object
    g = SingleGenotype(ens=ens,
                 ddg_dict=ddg_dict,
                 sites=[1],
                 mutations=["M1A"],
                 mutations_accumulated=["M1A"],
                 mut_energy=np.array([1,-1]))
    
    assert g.mutations_accumulated is g._mutations_accumulated
    assert np.array_equal(g.mutations_accumulated,["M1A"])


def test_SingleGenotype_sites(ens_test_data):
    
    ens = ens_test_data["ens"]
    ddg_dict = ens_test_data["ddg_dict"]

    # ---------------------------------
    # Create a genotype object
    g = SingleGenotype(ens=ens,
                 ddg_dict=ddg_dict,
                 sites=[1],
                 mutations=["M1A"],
                 mutations_accumulated=["M1A"],
                 mut_energy=np.array([1,-1]))
    
    assert g.sites is g._sites
    assert np.array_equal(g.sites,[1])



