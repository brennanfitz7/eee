
import eee
from eee.tools.plots import plt

import numpy as np

import json
import copy

def get_min_max(min_value,max_value,scalar=0.1,through_zero=False):

    # If min_value and max_value are both zero, expand them by scalar
    if min_value == max_value:
        if min_value == 0:
            offset = scalar/2
        else:
            offset = max_value*scalar

    else:
        offset = (max_value - min_value)*scalar

    min_value = min_value - offset
    max_value = max_value + offset
    
    # If we are missing zero, expand max or min appropriately
    if through_zero:
        if max_value < 0:
            max_value = 0 + (max_value - min_value)*scalar
        if min_value > 0:
            min_value = 0 - (max_value - min_value)*scalar

    return min_value, max_value



def traj_epistasis(sim_results,
                   key_traj,
                   ligand_name,
                   generate_plots=True):
    
    # Get the simulation
    so = sim_results["so"]

    # Genotypes seen (from most frequent to least)
    genotypes_seen = list(key_traj.keys())
    if "other" in genotypes_seen:
        genotypes_seen.remove("other")

    out = {}
    success = False
    for i in range(len(genotypes_seen)):

        try:
            genotype = sim_results["genotypes"].loc[int(genotypes_seen[i]),:]
        except IndexError:
            break

        if genotype.num_mutations == 2:

            # Get list of mutations
            mutations = genotype.mutations.split("/")

            # Case where we got two mutations at same site
            sites = [int(m[1:-1]) for m in mutations]
            if len(set(sites)) != 2:
                continue

            # Get mutations in correct order
            mutation_order = genotype.accum_mut.split("/")

            # Get some information about this genotype
            this_column = np.array(key_traj.loc[:,genotypes_seen[i]])
            max_freq = np.max(this_column)
            alive = this_column[-1] != 0
            sim_length = len(this_column)

            success = True

            break

    if not success:
        return None, None

    # Build initial output
    out = {}
    out["max_freq"] = max_freq
    out["alive"] = alive
    out["sim_length"] = sim_length

    sites = [int(m[1:-1]) for m in mutation_order]

    m1 = mutation_order[0]
    m2 = mutation_order[1]

    try:
        m1_ddg = so.gc.ddg_dict[sites[0]][m1]
        m2_ddg = so.gc.ddg_dict[sites[1]][m2]
    except KeyError:
        print(m1,m2)
        return None, None

    wt_fitness = so.fitness_from_energy(np.zeros(4))
    m1_fitness = so.fitness_from_energy(m1_ddg)
    m2_fitness = so.fitness_from_energy(m2_ddg)
    m12_fitness = so.fitness_from_energy(m1_ddg + m2_ddg)

    out["wt_fitness"] = wt_fitness
    out["m1_fitness"] = m1_fitness
    out["m2_fitness"] = m2_fitness
    out["m12_fitness"] = m12_fitness
    fit_ep = eee.analysis.epistasis.get_epistasis(wt_fitness,
                                                  m1_fitness,
                                                  m2_fitness,
                                                  m12_fitness)

    ligand_values = so.fc.ligand_dict[ligand_name]

    min_lig, max_lig = get_min_max(np.min(ligand_values),
                                   np.max(ligand_values))

    ligand_dict = copy.deepcopy(so.fc.ligand_dict)
    ligand_dict[ligand_name] = np.linspace(min_lig,max_lig,100)
    for lig in ligand_dict:
        if lig != ligand_name:
            ligand_dict[lig] = np.mean(ligand_dict[lig])
                                      
    m1_dict = {}
    m2_dict = {}
    m12_dict = {}
    for idx in range(len(m1_ddg)):
        m1_dict[so.ens.species[idx]] = m1_ddg[idx]
        m2_dict[so.ens.species[idx]] = m2_ddg[idx]
        m12_dict[so.ens.species[idx]] = m1_ddg[idx] + m2_ddg[idx]

    ens_ep = eee.analysis.epistasis.get_ensemble_epistasis(so.ens,
                                                           mut1_dict=m1_dict,
                                                           mut2_dict=m2_dict,
                                                           mut12_dict=m12_dict,
                                                           ligand_dict=ligand_dict)

    low_ligand = so.fc.ligand_dict[ligand_name][0]
    high_ligand = so.fc.ligand_dict[ligand_name][1]

    if generate_plots:

        fig, ax = plt.subplots(2,2,figsize=(12,12))
        fig.suptitle(f"{m1}/{m2}")
        ax[0,0].plot(ens_ep[ligand_name],ens_ep.fx_obs_00,lw=2,label="wt")
        ax[0,0].plot(ens_ep[ligand_name],ens_ep.fx_obs_10,lw=2,label=m1)
        ax[0,0].plot(ens_ep[ligand_name],ens_ep.fx_obs_01,lw=2,label=m2)
        ax[0,0].plot(ens_ep[ligand_name],ens_ep.fx_obs_11,lw=2,label=f"{m1}/{m2}")
        ax[0,0].legend()

        ax[0,0].plot((low_ligand,low_ligand),(0,1),'--',lw=1,color='gray')
        ax[0,0].plot((high_ligand,high_ligand),(0,1),'--',lw=1,color='gray')

        ax[0,0].set_ylabel("observable")

        ax[1,0].plot(ens_ep[ligand_name],ens_ep.fx_ep_mag,'-',lw=3,color='black')
        ax[1,0].plot(ens_ep[ligand_name],np.zeros(len(ens_ep[ligand_name])),'--',color='gray',zorder=-20)
        ax[1,0].plot((low_ligand,low_ligand),(-1,1),'--',lw=1,color='gray')
        ax[1,0].plot((high_ligand,high_ligand),(-1,1),'--',lw=1,color='gray')

        ax[1,0].set_ylabel("epistasis in fraction observable")
        ax[1,0].set_xlabel(f"{ligand_name} chemical potential")
                
        min_dG_ep, max_dG_ep = get_min_max(np.min(ens_ep.dG_ep_mag),
                                           np.max(ens_ep.dG_ep_mag),
                                           through_zero=True)

        ax[1,1].plot(ens_ep[ligand_name],ens_ep.dG_ep_mag,'-',lw=3,color='black')
        ax[1,1].plot(ens_ep[ligand_name],np.zeros(len(ens_ep[ligand_name])),'--',color='gray',zorder=-20)
        ax[1,1].plot((low_ligand,low_ligand),(min_dG_ep,max_dG_ep),'--',lw=1,color='gray')
        ax[1,1].plot((high_ligand,high_ligand),(min_dG_ep,max_dG_ep),'--',lw=1,color='gray') 

        ax[1,1].set_ylabel("epistasis in dG observable")
        ax[1,1].set_xlabel(f"{ligand_name} chemical potential")

        prop = dict(arrowstyle="-|>,head_width=0.4,head_length=0.8",
                    shrinkA=0,
                    shrinkB=0,
                    fc="black")

        path_taken = copy.deepcopy(prop)
        path_taken["lw"] = 2

        ax[0,1].annotate("", xy=(1,m1_fitness), xytext=(0,wt_fitness), arrowprops=path_taken)
        ax[0,1].annotate("", xy=(1,m2_fitness), xytext=(0,wt_fitness), arrowprops=prop)
        ax[0,1].annotate("", xy=(2,m12_fitness), xytext=(1,m1_fitness), arrowprops=path_taken)
        ax[0,1].annotate("", xy=(2,m12_fitness), xytext=(1,m2_fitness), arrowprops=prop)
        ax[0,1].set_xlim((-0.1,2.1))
        ax[0,1].set_ylim((-0.1,1.1))

        for i in range(2):
            for j in range(2):
                ax[i,j].spines['right'].set_visible(False)
                ax[i,j].spines['top'].set_visible(False)
        
        ax[0,1].spines["bottom"].set_visible(False)
        
        fig.tight_layout()

        plt.show()

    out["m1"] = m1
    out["m2"] = m2
    for k in m1_dict:
        out[f"m1_{k}_ddg"] = m1_dict[k]
    for k in m2_dict:
        out[f"m2_{k}_ddg"] = m2_dict[k]

    out["mag"] = fit_ep[0]
    out["sign1"] = fit_ep[1]
    out["sign2"] = fit_ep[2]
    out["ep_class"] = fit_ep[3]

    return out, ens_ep