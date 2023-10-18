
from eee.core.fitness.ff import ff_on
from eee.core.fitness.ff import ff_off
from eee.core.fitness.ff import ff_neutral
from eee.core.fitness.ff import ff_on_above
from eee.core.fitness.ff import ff_on_below

from eee.core.fitness.map_fitness_fcn import map_fitness_fcn

import pytest

def test_map_fitness_fcn(variable_types):
    
    assert map_fitness_fcn(value="on",return_as="string") == "on"
    assert map_fitness_fcn(value="on",return_as="function") == ff_on
    assert map_fitness_fcn(value=ff_on,return_as="string") == "on"
    assert map_fitness_fcn(value=ff_on,return_as="function") == ff_on

    assert map_fitness_fcn(value="off",return_as="string") == "off"
    assert map_fitness_fcn(value="off",return_as="function") == ff_off
    assert map_fitness_fcn(value=ff_off,return_as="string") == "off"
    assert map_fitness_fcn(value=ff_off,return_as="function") == ff_off

    assert map_fitness_fcn(value="neutral",return_as="string") == "neutral"
    assert map_fitness_fcn(value="neutral",return_as="function") == ff_neutral
    assert map_fitness_fcn(value=ff_neutral,return_as="string") == "neutral"
    assert map_fitness_fcn(value=ff_neutral,return_as="function") == ff_neutral

    assert map_fitness_fcn(value="on_above",return_as="string") == "on_above"
    assert map_fitness_fcn(value="on_above",return_as="function") == ff_on_above
    assert map_fitness_fcn(value=ff_on_above,return_as="string") == "on_above"
    assert map_fitness_fcn(value=ff_on_above,return_as="function") == ff_on_above

    assert map_fitness_fcn(value="on_below",return_as="string") == "on_below"
    assert map_fitness_fcn(value="on_below",return_as="function") == ff_on_below
    assert map_fitness_fcn(value=ff_on_below,return_as="string") == "on_below"
    assert map_fitness_fcn(value=ff_on_below,return_as="function") == ff_on_below

    for v in variable_types["not_str"]:
        if callable(v):
            continue
        
        if issubclass(type(v),type):
            continue

        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            map_fitness_fcn(v,return_as="string")    

    for v in variable_types["everything"]:
        if callable(v):
            continue

        if issubclass(type(v),type):
            continue

        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            map_fitness_fcn(v,return_as="function")    

    for v in variable_types["everything"]:
        if callable(v):
            continue

        if issubclass(type(v),type):
            continue

        print(v,type(v),flush=True)
        with pytest.raises(ValueError):
            map_fitness_fcn("on",return_as=v)    
