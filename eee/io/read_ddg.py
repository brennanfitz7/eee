"""
Read a ddg spreadsheet file.
"""
import numpy as np

from eee.io.read_dataframe import read_dataframe

def read_ddg(ddg_df):
    """
    Load a ddg file, removing any self-to-self mutations. 

    Parameters
    ----------
    ddg_df : str OR pandas.DataFrame
        spreadsheet file OR loaded dataframe. Should have a "mut" column 
        (formatted like A21A, Q45L, etc.) and columns for each species in the
        ensemble. The values in the species columns are the predicted energetic
        effect of that mutation on that ensemble species. 
        
    Returns
    -------
    df : pandas.DataFrame
        pandas dataframe with columns holding 'mut', 'site' (i.e., the 21 in 
        A21P), and then columns for the predicted ddG for each species. Any 
        columns except "mut" and "site" will be untouched by the function. 
        
    Example
    -------
    Here is an example ddg_file for an ensemble with the species erg, erg-oht, 
    and erg-oht-pep.
    
    ..code-block::
    
        mut,erg,erg-oht,erg-oht-pep
        Y1A,-0.1510000000000673,-0.12999999999999545,0.4900000000000091
        Y1C,0.4880000000000564,0.17599999999993088,2.2889999999999873
        Y1D,-2.433999999999969,-2.3319999999999936,-1.6990000000000691
        Y1F,-4.645999999999958,-0.4049999999999727,3.1989999999999554
        
    
    """
    
    # Read dataframe (function handles xlsx, csv, tsv, or pandas.DataFrame)
    df = read_dataframe(ddg_df)

    # Check to make sure there is a mut column
    if "mut" not in df.columns:
        err = f"ddg_df does not have a 'mut' column.\n"
        raise ValueError(err)

    # Find the wildtype entries (i.e., Q45Q). 
    wt_mask = np.array([m[0] == m[-1] for m in df.mut],dtype=bool)
    df = df.loc[np.logical_not(wt_mask),:]

    # Extract sites
    df["site"] = [int(m[1:-1]) for m in list(df.mut)]

    # Make sure the dataframe starts with site, mut, then other stuff
    columns = list(df.columns)
    columns.remove("site")
    columns.remove("mut")
    columns.insert(0,"site")
    columns.insert(1,"mut")

    # Reorder columns
    df = df.loc[:,columns]
    
    return df


