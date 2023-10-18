"""
Functions to analyze epistasis.
"""

import numpy as np
import pandas as pd
from tqdm.auto import tqdm


def get_epistasis(m00,m10,m01,m11):
    """
    Get epistasis for values given in m00, m10, m01, and m11. 
    
    Parameters
    ----------
    m00 : float or np.array
        value of observable without mutations. can be a float or an array of
        floats. If an array, it must have the same length as arrays for m10, 
        m01, and m11.
    m10 : float or np.array
        value of observable with mutation 1
    m01 : float or np.array
        value of observable with mutation 2
    m11 : float o rnp.array
        value of observable with both mutations 1 and 2
        
    Returns
    -------
    mag : float or np.array
        magnitude of the epistasis
    sign1 : bool or np.array
        whether mutation 1 changes sign due to mutation 2
    sign2 : bool or np.array
        whether mutation 2 changes sign due to mutation 1
    ep_class : str or np.array
        type of epistasis seen. None (no epistasis), "mag" (magnitude), 
        "sign" (sign), or "recip" (reciprocal sign).
    """

    # magnitude of epistasis (signed)
    mag = (m11 - m10) - (m01 - m00)
    
    # Sign of mutation 1 (will be False if mutation effect has same
    # sign in both backgrounds; True if opposite signs)
    sign1 = (m11 - m01)/(m10 - m00)
    sign1 = sign1 < 0

    # Sign of mutation 1 (will be False if mutation effect has same
    # sign in both backgrounds; True if opposite signs)
    sign2 = (m11 - m10)/(m01 - m00)
    sign2 = sign2 < 0

    # Is epistasis > 0?
    is_epistasis = np.logical_not(np.isclose(mag,0))
    
    # Some kind of sign epistasis
    is_sign = np.logical_or(sign1,sign2)

    # Magnitude if not sign
    is_mag = np.logical_not(is_sign)
    
    # Separate reciprocal from simple sign epistasis
    is_recip = np.logical_and(sign1,sign2)
    
    is_sign = np.logical_and(is_sign,np.logical_not(is_recip))

    # Filter all classes of epistasis based on magnitude
    is_mag = np.logical_and(is_mag,is_epistasis)
    is_sign = np.logical_and(is_sign,is_epistasis)
    is_recip = np.logical_and(is_recip,is_epistasis)
    
    # Record classes. If iterable, create numpy array. If single value, return
    # single values
    if hasattr(mag,'__iter__'):
        ep_class = np.array([None for _ in range(len(mag))])
        ep_class[is_mag] = "mag"
        ep_class[is_sign] = "sign"
        ep_class[is_recip] = "recip"
    else:
        if is_recip:
            ep_class = "recip"
        elif is_sign:
            ep_class = "sign"
        elif is_mag:
            ep_class = "mag"
        else:
            ep_class = None

    return mag, sign1, sign2, ep_class


def get_ensemble_epistasis(ens,
                           mut1_dict=None,
                           mut2_dict=None,
                           mut12_dict=None,
                           ligand_dict=None,
                           temperature=298.15):
    """
    Get the epistasis between two mutations across different chemical 
    potentials.
    
    Parameters
    ----------
    ens : eee.core.Ensemble 
        ensemble instance whose species whose names match the keys in the
        mut_dicts and whose ligands match the keys in ligand_dict
    mut_dict1 : dict
        dictionary holding effects of mutation 1 on different species. 
    mut_dict2 : dict
        dictionary holding effects of mutation 2 on different species. 
    mut_dict12 : dict
        dictionary holding combined effects of mutations 1 and 2 on
        different species. 
    ligand_dict : dict
        dictionary of chemical potentials
    temperature : float, default=298.15
        temperature in Kelvin

    Returns
    -------
    df : pandas.DataFrame
        pandas dataframe with fx_obs for each genotype (00,10,01,11), 
        the magnitude, sign, and class (mag, sign, reciprocal 
        sign) of epistasis in fx_obs, dG_obs for each genotype, and then
        the epistasis in dG_obs. These are reported as a function of the 
        species concentrations in ligand_dict. 
    """

    # Calculate observables for each genotype
    df_00 = ens.get_obs(mut_energy=None,
                        ligand_dict=ligand_dict,
                        temperature=temperature)

    df_10 = ens.get_obs(mut_energy=mut1_dict,
                        ligand_dict=ligand_dict,
                        temperature=temperature)

    df_01 = ens.get_obs(mut_energy=mut2_dict,
                        ligand_dict=ligand_dict,
                        temperature=temperature)

    df_11 = ens.get_obs(mut_energy=mut12_dict,
                        ligand_dict=ligand_dict,
                        temperature=temperature)

    # Create dataframe
    columns = ens.ligands[:]
    columns.insert(0,"temperature")
    df = df_00.loc[:,columns]

    # Epistasis in fx_obs
    df["fx_obs_00"] = df_00.loc[:,"fx_obs"]
    df["fx_obs_10"] = df_10.loc[:,"fx_obs"]
    df["fx_obs_01"] = df_01.loc[:,"fx_obs"]
    df["fx_obs_11"] = df_11.loc[:,"fx_obs"]

    ep_mag, ep_sign1, ep_sign2, ep_class = get_epistasis(df_00.loc[:,"fx_obs"],
                                                         df_10.loc[:,"fx_obs"],
                                                         df_01.loc[:,"fx_obs"],
                                                         df_11.loc[:,"fx_obs"])
        
    df["fx_ep_mag"] = ep_mag
    df["fx_ep_sign1"] = ep_sign1
    df["fx_ep_sign2"] = ep_sign2
    df["fx_ep_class"] = ep_class

    # Epistasis in dG_obs
    df["dG_obs_00"] = df_00.loc[:,"dG_obs"]
    df["dG_obs_10"] = df_10.loc[:,"dG_obs"]
    df["dG_obs_01"] = df_01.loc[:,"dG_obs"]
    df["dG_obs_11"] = df_11.loc[:,"dG_obs"]

    ep_mag, ep_sign1, ep_sign2, ep_class = get_epistasis(df_00.loc[:,"dG_obs"],
                                                         df_10.loc[:,"dG_obs"],
                                                         df_01.loc[:,"dG_obs"],
                                                         df_11.loc[:,"dG_obs"])
        
    df["dG_ep_mag"] = ep_mag
    df["dG_ep_sign1"] = ep_sign1
    df["dG_ep_sign2"] = ep_sign2
    df["dG_ep_class"] = ep_class

    return df


def get_all_pairs_epistasis(ens,ddg_df,ligand_dict,get_only=None):
    """
    Given an ensemble and a mutation ddG file, calculate epistasis for all 
    pairs of mutations.
    
    Parameters
    ----------
    ens : Ensemble instance
        ensemble instance holding species whose names match the column names
        in a ddg csv file
    ddg_file : pandas.DataFrame
        dataframe holding the energetic effects mutations on the species
        in the ensemble. 
    ligand_dict : dict
        dictionary holding chemical potentials at which to make the calculation.
        See Ensemble docs for details.
    get_only : int, optional
        if set, take only get_only lines from the ddg file.
        
    Returns
    -------
    out_df : pandas.DataFrame
        dataframe with six columns: m1 (mutation 1), m2 (mutation 2), s1 (site 1),
        s2 (site 2), ep_mag (magnitude of the observed epistasis), ep_class 
        (class of epistasis: mag, sign, recip_sign). 
    """


    if get_only is not None:
        if get_only < len(ddg_df.index):
            idx = np.random.choice(ddg_df.index,size=get_only,replace=False)
            ddg_df = ddg_df.loc[idx,:]
    
    # Get names of molecular species
    species = ens.species 

    # Create output dictionary (source for dataframe)
    columns = ["m1","m2","s1","s2","ep_mag","ep_class"]
    out = dict([(c,[]) for c in columns])
    
    # Create status bar
    N = len(ddg_df.index)
    total = N*(N-1)//2
    with tqdm(total=total) as p:

        # Go over all entries
        for i in range(len(ddg_df.index)):
            
            # i is wildtype
            if np.sum(ddg_df.loc[ddg_df.index[i],"is_wt"]) > 0:
                p.update(1)
                continue

            # Mutational effects at site i
            mut1_dict = {}
            for s in species:
                mut1_dict[s] = ddg_df.loc[ddg_df.index[i],s]

            # Paired entries
            for j in range(i+1,len(ddg_df.index)):

                # i and j are at same site
                if ddg_df.loc[ddg_df.index[i],"site"] == ddg_df.loc[ddg_df.index[j],"site"]:
                    p.update(1)
                    continue
                
                # j is wildtype
                if np.sum(ddg_df.loc[ddg_df.index[j],"is_wt"]) > 0:
                    p.update(1)
                    continue

                # Mutational effects at site j
                mut2_dict = {}
                for s in species:
                    mut2_dict[s] = ddg_df.loc[ddg_df.index[j],s]

                # Get epistasis
                df_ep = ens.get_epistasis(mut1_dict=mut1_dict,
                                          mut2_dict=mut2_dict,
                                          ligand_dict=ligand_dict)

                # Get mag and class
                ep_mag = np.max(df_ep.dG_ep_mag)
                ep_class = df_ep["dG_ep_class"].iloc[np.argmax(df_ep.dG_ep_mag)]

                # Record output
                out["m1"].append(ddg_df.loc[ddg_df.index[i],"mut"])
                out["m2"].append(ddg_df.loc[ddg_df.index[j],"mut"])
                out["s1"].append(ddg_df.loc[ddg_df.index[i],"site"])
                out["s2"].append(ddg_df.loc[ddg_df.index[j],"site"])
                out["ep_mag"].append(ep_mag)
                out["ep_class"].append(ep_class)
                
                p.update(1)

    out_df = pd.DataFrame(out)

    return out_df

def summarize_epistasis(df,cutoff=None):
    """
    Summarize the epistasis in a site-by-site manner.
    
    Parameters
    ----------
    df : pandas.DataFrame
        output from get_all_pairs_epistasis
    cutoff : float, optional
        only keep epistasis whose absolute mangitude is above cutoff
    
    Returns
    -------
    summary_df : pandas.DataFrame
        dataframe with columns holding summary information for each site. The 
        columns are the mean absolute value epistatic magnitude, the mean epistatic
        magntitutde, the standard deviation of the epistaitc magnitude, the 
        number of times the site has no epistaiss, the number of time the site has
        magnitude epistasis, the number of times the site has reciprocal epistasis,
        and the number of times the site has sign epistasis. 
    """

    columns = ["site",
               "abs_mag","avg_mag","std_mag",
               "ep_none","ep_mag","ep_recip","ep_sign"]
    out = dict([(c,[]) for c in columns])

    s1_sites = np.unique(df.s1)
    s2_sites = np.unique(df.s2)
    all_sites = list(s1_sites)
    all_sites.extend(s2_sites)
    all_sites = list(set(all_sites))
    all_sites = np.array(all_sites)

    for s in all_sites:

        site_mask = np.logical_or(df.s1 == s,df.s2 == s)
        site_df = df.loc[site_mask,:]
        
        # Drop epistasis below the cutoff
        if cutoff is not None:
            above_cutoff = np.abs(site_df.ep_mag) > cutoff
            site_df = site_df.loc[above_cutoff,:]

        abs_mag = np.mean(np.abs(site_df.ep_mag))
        avg_mag = np.mean(site_df.ep_mag)
        std_mag = np.std(site_df.ep_mag)
        nothing_mask = pd.isnull(df.ep_class)

        classy_sites = site_df.loc[np.logical_not(nothing_mask),"ep_class"]
        class_bins, class_counts = np.unique(classy_sites,return_counts=True)

        count_dict = dict(zip(class_bins,class_counts))
        for c in ["mag","recip","sign"]:
            if c not in count_dict:
                count_dict[c] = 0

        out["site"].append(s)
        out["abs_mag"].append(abs_mag)
        out["avg_mag"].append(avg_mag)
        out["std_mag"].append(std_mag)
        out["ep_none"].append(np.sum(nothing_mask))
        out["ep_mag"].append(count_dict["mag"])
        out["ep_recip"].append(count_dict["recip"])
        out["ep_sign"].append(count_dict["sign"])


    summary_df = pd.DataFrame(out)

    return summary_df


