import pytest

from eee.analysis.cluster_ddg import _clean_up_axes
from eee.analysis.cluster_ddg import cluster_ddg

import numpy as np
import pandas as pd

import os
import glob

def test__clean_up_axes():
    
    # Basic calc
    out = _clean_up_axes(min_value=0.0,
                         max_value=1.0,
                         num_ticks=5,
                         scalar=0)
    
    assert np.isclose(out[0],0)
    assert np.isclose(out[1],1)
    assert np.array_equal(out[2],np.arange(5)*0.25)

    # Check span change
    out = _clean_up_axes(min_value=-1.0,
                         max_value=1.0,
                         num_ticks=5,
                         scalar=0)
    
    assert np.isclose(out[0],-1)
    assert np.isclose(out[1],1)
    assert np.array_equal(out[2],np.arange(-2,3)*0.5)

    # Check min/max flip
    out = _clean_up_axes(min_value=1.0,
                         max_value=-1.0,
                         num_ticks=5,
                         scalar=0)
    
    assert np.isclose(out[0],-1)
    assert np.isclose(out[1],1)
    assert np.array_equal(out[2],np.arange(-2,3)*0.5)

    with pytest.raises(ValueError):
        out = _clean_up_axes(min_value=0.0,
                             max_value=0.0,
                             num_ticks=5,
                             scalar=0)

    # test expansion
    # Basic calc
    out = _clean_up_axes(min_value=0.0,
                         max_value=1.0,
                         num_ticks=5,
                         scalar=0.05)
    
    assert np.isclose(out[0],-0.05)
    assert np.isclose(out[1],1.05)


    out = _clean_up_axes(min_value=-1.0,
                         max_value=1.0,
                         num_ticks=5,
                         scalar=0.05)
    
    assert np.isclose(out[0],-1.1)
    assert np.isclose(out[1],1.1)


    out = _clean_up_axes(min_value=-2.0,
                         max_value=-1.0,
                         num_ticks=5,
                         scalar=0.05)
    
    assert np.isclose(out[0],-2.05)
    assert np.isclose(out[1],-0.95)

def test_cluster_ddg(tmpdir,test_ddg):
    
    current_dir = os.getcwd()
    os.chdir(tmpdir)

    should_exist = ["yo_kmeans-elbow-plot.pdf",
                    "yo_clusters.pdf"]
    for e in should_exist:
        assert not os.path.isfile(e)

    out = cluster_ddg(test_ddg["ddg.csv"],
                      max_num_clusters=10,
                      elbow_cutoff=0.1,
                      write_prefix="yo",
                      exclude_columns=None)
    
    for e in should_exist:
        assert os.path.isfile(e)
    
    assert len(np.unique(out["cluster"])) == 4
    for g in glob.glob("yo*"):
        os.remove(g)
    
    out = cluster_ddg(test_ddg["ddg.csv"],
                      max_num_clusters=1,
                      elbow_cutoff=0.1,
                      write_prefix="yo",
                      exclude_columns=None)
    
    assert len(np.unique(out["cluster"])) == 1
    for g in glob.glob("yo*"):
        os.remove(g)

    out = cluster_ddg(test_ddg["ddg.csv"],
                      max_num_clusters=2,
                      elbow_cutoff=0.1,
                      write_prefix="yo",
                      exclude_columns=None)
    
    assert len(np.unique(out["cluster"])) == 2
    for g in glob.glob("yo*"):
        os.remove(g)
    
    with pytest.raises(ValueError):
        out = cluster_ddg(test_ddg["ddg.csv"],
                        max_num_clusters=0,
                        elbow_cutoff=0.1,
                        write_prefix="yo",
                        exclude_columns=None)

    # change elbow cutoff
    out = cluster_ddg(test_ddg["ddg.csv"],
                      max_num_clusters=10,
                      elbow_cutoff=0.0000001,
                      write_prefix="yo",
                      exclude_columns=None)
    
    assert len(np.unique(out["cluster"])) == 10
    for g in glob.glob("yo*"):
        os.remove(g)


    # Change write prefix
    should_exist = ["yo2_kmeans-elbow-plot.pdf",
                    "yo2_clusters.pdf"]
    for e in should_exist:
        assert not os.path.isfile(e)

    out = cluster_ddg(test_ddg["ddg.csv"],
                      max_num_clusters=10,
                      elbow_cutoff=0.1,
                      write_prefix="yo2",
                      exclude_columns=None)
    
    for e in should_exist:
        assert os.path.isfile(e)
    
    assert len(np.unique(out["cluster"])) == 4
    for g in glob.glob("yo*"):
        os.remove(g)
    
    # Change write prefix to None
    out = cluster_ddg(test_ddg["ddg.csv"],
                      max_num_clusters=10,
                      elbow_cutoff=0.1,
                      write_prefix=None,
                      exclude_columns=None)
    assert len(glob.glob("*.pdf")) == 0
    assert len(np.unique(out["cluster"])) == 4
    
    # Test exclude_columns
    df = pd.read_csv(test_ddg["ddg.csv"])
    df["stupid"]  = 1
    exclude_columns = ["stupid"]

    out = cluster_ddg(test_ddg["ddg.csv"],
                      max_num_clusters=10,
                      elbow_cutoff=0.1,
                      write_prefix="yo2",
                      exclude_columns=exclude_columns)
    
    assert "stupid" not in out.columns
    for g in glob.glob("yo*"):
        os.remove(g)
    

    os.chdir(current_dir)