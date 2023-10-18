"""
Write a tree in newick format.
"""

from eee._private.check.standard import check_bool
from eee._private.check.standard import check_int

import ete3

import os

def write_tree(T,
               out_file=None,
               fmt=3,
               overwrite=False):
    """
    Write out an ete3.Tree as a newick format. 

    Parameters
    ----------
    T : ete3.TreeNode
        ete3 tree
    out_file : str, optional
        output file. If defined, write the newick string the file.
    fmt : int, default=3
        tree writing format. See the docstring for :code:`read_tree` for details.
        default writes node names, ancestor names, and branch lengths. 
    overwrite : bool, default=False
        whether or not to overwrite an existing file

    Returns
    -------
    tree : str
        Newick string representation of the output tree(s)
    """
    
    # --------------------------------------------------------------------------
    # Parameter sanity checking

    if not issubclass(type(T),ete3.TreeNode):
        err = "\nT must be an ete3.Tree instance\n\n"
        raise ValueError(err)
    
    fmt = check_int(fmt,
                    variable_name="fmt",
                    minimum_allowed=0,
                    maximum_allowed=100)
    if fmt > 9 and fmt != 100:
        err = "\nfmt must be a value between 0 and 9, or 100.\n\n"
        raise ValueError(err)
    
    if out_file is not None:

        if not issubclass(type(out_file),str):
            err = "\nout_file must be a string pointing to a file to write out\n\n"
            raise ValueError(err)
        
        overwrite = check_bool(overwrite,"overwrite")

        if os.path.exists(out_file):
            if os.path.isfile(out_file):
                if overwrite:
                    os.remove(out_file)
                else:
                    err = f"\nout_file '{out_file}' exists. Either delete or set overwrite to True\n\n"
                    raise FileExistsError(err)
            else:
                err = f"\nout_file '{out_file}' exists but is a directory. Cannot write output.\n\n"
                raise FileExistsError(err)

    as_string = T.write(format=fmt)

    if out_file is not None:
        with open(out_file,'w') as f:
            f.write(as_string)
    
    return as_string