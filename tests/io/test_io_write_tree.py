import pytest

from eee.io.read_tree import read_tree
from eee.io.write_tree import write_tree

import os

def test_write_tree(newick_files,variable_types,tmpdir):

    cwd = os.getcwd()
    os.chdir(tmpdir)

    T = read_tree(newick_files["simple.newick"])

    formats = list(range(10))
    formats.append(100)
    for f in formats:
        T_to_write = T.copy()
        T2 = write_tree(T_to_write,
                         out_file="test.newick",
                         fmt=f)
        assert os.path.exists("test.newick")
        os.remove("test.newick")

    # make sure it is sending out tree
    T_out = write_tree(T,fmt=3,out_file=None)
    assert issubclass(type(T_out),str)
    assert len(T_out) == 50

    T_short_out = write_tree(T,fmt=100,out_file=None)
    assert issubclass(type(T_short_out),str)

    # Only topology tree will be a shorter string than the tree with branch 
    # lengths and labels
    assert len(T_short_out) < len(T_out)

    # Try to write to existing
    with open('test.newick','w') as f:
        f.write("block")

    with pytest.raises(FileExistsError):
        newick = write_tree(T,out_file="test.newick")
                
    # Overwrite
    _ = write_tree(T,out_file="test.newick",overwrite=True)
    os.remove("test.newick")

    # Try to write onto a directory
    os.mkdir("test.newick")
    with pytest.raises(FileExistsError):
        _ = write_tree(T,out_file="test.newick",overwrite=True)
    
    for v in variable_types["everything"]:
        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            write_tree(T=v)

    for v in variable_types["everything"]:
        if issubclass(type(v),str):
            continue
        if v is None:
            continue
        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            write_tree(T=T,out_file=v)

    for v in variable_types["not_ints_or_coercable"]:
        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            write_tree(T=T,fmt=v,out_file=None) 
    
    for v in [-1,10,99,101]:
        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            write_tree(T=T,fmt=v,out_file=None) 

    os.chdir(cwd)