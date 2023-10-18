
import pytest

from eee._private.check.compare_dict import compare_dict

import numpy as np
import pandas as pd

def test_compare_dict(variable_types):

    assert compare_dict({},{}) == True
    assert compare_dict("test","test") == False
    assert compare_dict({"test":1},"test") == False
    assert compare_dict("test",{"test":1}) == False
    assert compare_dict({"test":1},{"test":1}) == True
    assert compare_dict({"test":1},{"test":2}) == False
    assert compare_dict({"test":1},{"taste":1}) == False
    assert compare_dict({"test":1,"this":2},{"taste":1,"this":2}) == False

    # Identical complicated dict
    x1 = {"test":{"this":1,"out":2},
        "yo":{"a":3,"b":4}}

    x2 = {"test":{"this":1,"out":2},
        "yo":{"a":3,"b":4}}
    assert compare_dict(x1,x2) == True

    # Different integer
    x1 = {"test":{"this":1,"out":2},
        "yo":{"a":3,"b":4}}

    x2 = {"test":{"this":1,"out":2},
        "yo":{"a":3,"b":5}}
    assert compare_dict(x1,x2) == False

    # Identical set and arrays in complicated dict
    x1 = {"test":{"this":np.array([1,2,3]),"out":2},
        "yo":{"a":set([1,2,3]),"b":4}}

    x2 = {"test":{"this":np.array([1,2,3]),"out":2},
        "yo":{"a":set([1,2,3]),"b":4}}
    assert compare_dict(x1,x2) == True

    # Different set in complicated dict
    x1 = {"test":{"this":np.array([1,2,3]),"out":2},
        "yo":{"a":set([1,2,3]),"b":4}}

    x2 = {"test":{"this":np.array([1,2,3]),"out":2},
        "yo":{"a":set([1,2,3,4]),"b":4}}
    assert compare_dict(x1,x2) == False

    # Set versus list in complicated dict
    x1 = {"test":{"this":np.array([1,2,3]),"out":2},
        "yo":{"a":[1,2,3],"b":4}}

    x2 = {"test":{"this":np.array([1,2,3]),"out":2},
        "yo":{"a":set([1,2,3]),"b":4}}
    assert compare_dict(x1,x2) == False

    # Array versus list (identical entries) in complicated dict
    x1 = {"test":{"this":np.array([1,2,3]),"out":2},
        "yo":{"a":set([1,2,3]),"b":4}}

    x2 = {"test":{"this":[1,2,3],"out":2},
        "yo":{"a":set([1,2,3]),"b":4}}
    assert compare_dict(x1,x2) == True

    # Array versus list (different entries) in complicated dict
    x1 = {"test":{"this":np.array([1,2,4]),"out":2},
        "yo":{"a":set([1,2,3]),"b":4}}

    x2 = {"test":{"this":[1,2,3],"out":2},
        "yo":{"a":set([1,2,3]),"b":4}}
    assert compare_dict(x1,x2) == False

    x1 = {"test":pd.DataFrame({"x":[1,2,3]})}
    x2 = {"test":pd.DataFrame({"x":[1,2,3]})}
    assert compare_dict(x1,x2) == True

    x1 = {"test":pd.DataFrame({"x":[1,2,3]})}
    x2 = {"test":pd.DataFrame({"y":[1,2,3]})}
    assert compare_dict(x1,x2) == False

    x1 = {"test":pd.DataFrame({"x":[1,2]})}
    x2 = {"test":pd.DataFrame({"x":[1,2,3]})}
    assert compare_dict(x1,x2) == False

    x1 = {"test":pd.DataFrame({"x":[1,2,3],"y":[1,3,4]})}
    x2 = {"test":pd.DataFrame({"x":[1,2,3]})}
    assert compare_dict(x1,x2) == False

    for v in variable_types["not_dict"]:
        print(v,type(v),flush=True)
        assert compare_dict(v,v) == False

    assert compare_dict({"test":np.nan},{"test":np.nan}) == False
    assert compare_dict({"test":pd.NA},{"test":pd.NA}) == False

    for v in variable_types["everything"]:

        print(v,type(v),flush=True)

        # hacked test. check for nan and null values for everything except 
        # iterables. 
        if not hasattr(v,"__iter__"):

            try:
                if np.isnan(v):
                    continue
            except (TypeError,ValueError):
                pass

            try:
                if pd.isnull(v):
                    continue
            except (TypeError,ValueError):    
                pass

        assert compare_dict({"test":v},{"test":v}) == True
        
    # A few more case checks
    x1 = {"test":{"this":1}}
    x2 = {"test":"this"}
    assert compare_dict(x1,x2) == False

    x1 = {"test":set(["this"])}
    x2 = {"test":["this"]}
    assert compare_dict(x1,x2) == False

    x1 = {"test":pd.DataFrame({"this":[1]})}
    x2 = {"test":{"this":1}}
    assert compare_dict(x1,x2) == False

    x1 = {"test":pd.DataFrame({"this":[1,2]})}
    x2 = {"test":pd.DataFrame({"this":[1,2,3]})}
    assert compare_dict(x1,x2) == False

    x1 = {"test":[1]}
    x2 = {"test":1}
    assert compare_dict(x1,x2) == False

    x1 = {"test":[1]}
    x2 = {"test":[1,2]}
    assert compare_dict(x1,x2) == False

    x1 = {"test":np.nan}
    x2 = {"test":1}
    assert compare_dict(x1,x2) == False

    x1 = {"test":1}
    x2 = {"test":np.nan}
    assert compare_dict(x1,x2) == False

    x1 = {"test":np.nan}
    x2 = {"test":np.nan}
    assert compare_dict(x1,x2) == False

    x1 = {"test":pd.NA}
    x2 = {"test":1}
    assert compare_dict(x1,x2) == False

    x1 = {"test":1}
    x2 = {"test":pd.NA}
    assert compare_dict(x1,x2) == False

    x1 = {"test":pd.NA}
    x2 = {"test":pd.NA}
    assert compare_dict(x1,x2) == False