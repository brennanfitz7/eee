"""
Flexible function for reading a dataframe. Handles xlsx, csv, tsv, and DataFrame
inputs.
"""

import pandas as pd
import numpy as np

def read_dataframe(df,remove_extra_index=True):
    """
    Read a spreadsheet. Handles .csv, .tsv, .xlsx/.xls. If extension is
    not one of these, attempts to parse text as a spreadsheet using
    :code:`pandas.read_csv(sep=None)`.

    Parameters
    ----------
    input : pandas.DataFrame or str
        either a pandas dataframe OR the filename to read in.
    remove_extra_index : bool, default=True
        look for any 'Unnamed: xx' columns that pandas writes out for
        pandas.to_csv(index=True) and, if found, drop.

    Returns
    -------
    pandas.DataFrame
        validated dataframe
    """

    # If this is a string, try to load it as a file
    if issubclass(type(df),str):

        filename = df

        ext = filename.split(".")[-1].strip().lower()

        if ext in ["xlsx","xls"]:
            df = pd.read_excel(filename)
        elif ext == "csv":
            df = pd.read_csv(filename,sep=",")
        elif ext == "tsv":
            df = pd.read_csv(filename,sep="\t")
        else:
            # Fall back -- try to guess delimiter
            df = pd.read_csv(filename,sep=None,engine="python")

    # If this is a pandas dataframe, work in a copy of it.
    elif issubclass(type(df),pd.DataFrame):
        df = df.copy()

    # Otherwise, fail
    else:
        err = f"\n\n'df' {df} not recognized. Should be the filename of\n"
        err += "spreadsheet or a pandas dataframe.\n"
        raise ValueError(err)

    # Look for extra index columns that pandas writes out (in case user wrote out
    # pandas frame manually, then re-read). Looks for columns that start with 
    # Unnamed and have data type int. 
    if remove_extra_index:

        to_drop = []

        possible_extra = [c for c in df.columns if c.startswith("Unnamed:")]
        for column in possible_extra:
            possible_index = df.loc[:,column]
            if np.issubdtype(possible_index.dtypes,int):
                to_drop.append(column)
        
        if len(to_drop) > 0:
            df = df.drop(columns=to_drop)

    return df