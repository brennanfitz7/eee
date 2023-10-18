from eee.analysis.epistasis import get_epistasis
from eee.analysis.epistasis import get_ensemble_epistasis
from eee.analysis.epistasis import get_all_pairs_epistasis
from eee.analysis.epistasis import summarize_epistasis

import numpy as np

def test_get_epistasis():
    
    # No epistasis
    mag, sign1, sign2, ep_class = get_epistasis(m00=0,
                                                m10=1,
                                                m01=1,
                                                m11=2)
    assert mag == 0
    assert sign1 == False
    assert sign2 == False
    assert ep_class is None

    # Magnitude
    mag, sign1, sign2, ep_class = get_epistasis(m00=0,
                                                m10=1,
                                                m01=1,
                                                m11=3)
    assert mag == 1
    assert sign1 == False
    assert sign2 == False
    assert ep_class == "mag"

    mag, sign1, sign2, ep_class = get_epistasis(m00=0,
                                                m10=1,
                                                m01=1,
                                                m11=1)
    assert mag == -1
    assert sign1 == False
    assert sign2 == False
    assert ep_class == "mag"
    
    # Sign
    mag, sign1, sign2, ep_class = get_epistasis(m00=0,
                                                m10=-1,
                                                m01=1,
                                                m11=2)
    assert mag == 2
    assert sign1 == True
    assert sign2 == False
    assert ep_class == "sign"

    # Sign
    mag, sign1, sign2, ep_class = get_epistasis(m00=0,
                                                m10=1,
                                                m01=-1,
                                                m11=2)
    assert mag == 2
    assert sign1 == False
    assert sign2 == True
    assert ep_class == "sign"

    # Reciprocal sign
    mag, sign1, sign2, ep_class = get_epistasis(m00=0,
                                                m10=-1,
                                                m01=-1,
                                                m11=1)
    assert mag == 3
    assert sign1 == True
    assert sign2 == True
    assert ep_class == "recip"


    # Iterable with all possible
    m00 = np.array([0,0,0,0,0,0])
    m10 = np.array([1,1,1,-1,1,-1])
    m01 = np.array([1,1,1,1,-1,-1])
    m11 = np.array([2,3,1,2,2,1])

    mag, sign1, sign2, ep_class = get_epistasis(m00=m00,
                                                m10=m10,
                                                m01=m01,
                                                m11=m11)

    assert np.array_equal(mag,[0,1,-1,2,2,3])
    assert np.array_equal(sign1,[0,0,0,1,0,1])
    assert np.array_equal(sign2,[0,0,0,0,1,1])
    assert np.array_equal(ep_class,[None,"mag","mag","sign","sign","recip"])



def test_get_ensemble_epistasis():
    pass

def test_get_all_pairs_epistasis():
    pass

def test_summarize_epistasis():
    pass