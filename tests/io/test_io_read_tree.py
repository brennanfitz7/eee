import pytest

from eee.io.read_tree import read_tree

import ete3

import os

def test_read_tree(newick_files,variable_types,tmpdir):

    cwd = os.getcwd()
    os.chdir(tmpdir)

    T = read_tree(newick_files["simple.newick"])
    assert issubclass(type(T),ete3.TreeNode)

    T2 = read_tree(T)
    assert T2 is T

    for v in variable_types["everything"]:
        if issubclass(type(v),str):
            continue
        if v is None:
            continue
        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            read_tree(tree=v)

    for v in variable_types["not_ints_or_coercable"]:
        if v is None:
            continue
        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            read_tree(tree=newick_files["simple.newick"],
                      fmt=v) 
    
    for v in [-1,10,99,101]:
        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            read_tree(tree=newick_files["simple.newick"],
                      fmt=v) 


    os.chdir(cwd)
    