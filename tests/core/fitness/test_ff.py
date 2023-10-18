from eee.core.fitness.ff import ff_on
from eee.core.fitness.ff import ff_off
from eee.core.fitness.ff import ff_neutral
from eee.core.fitness.ff import ff_on_above
from eee.core.fitness.ff import ff_on_below
from eee.core.fitness.ff import FF_AVAILABLE

def test_ff_on():
    assert ff_on(1) == 1
    assert ff_on(0) == 0

    # test robustness to extra args
    assert ff_on(1,5,some_random_kwarg="test") == 1
    assert ff_on(0,5,some_random_kwarg="test") == 0

def test_ff_off():
    assert ff_off(1) == 0
    assert ff_off(0) == 1

    # test robustness to extra args
    assert ff_off(1,5,some_random_kwarg="test") == 0
    assert ff_off(0,5,some_random_kwarg="test") == 1

def test_ff_neutral():
    assert ff_neutral(1) == 1
    assert ff_neutral(0) == 1

    # test robustness to extra args
    assert ff_neutral(0,5,some_random_kwarg="test") == 1
    assert ff_neutral(1,5,some_random_kwarg="test") == 1

def test_ff_on_above():
    assert ff_on_above(5,threshold=4) == 1
    assert ff_on_above(5,threshold=5) == 1
    assert ff_on_above(3,threshold=4) == 0

    # test robustness to extra args
    assert ff_on_above(5,5,threshold=4,some_random_kwarg="test") == 1
    assert ff_on_above(3,5,threshold=4,some_random_kwarg="test") == 0

def test_ff_on_below():
    assert ff_on_below(5,threshold=4) == 0
    assert ff_on_below(5,threshold=5) == 1
    assert ff_on_below(3,threshold=4) == 1

    # test robustness to extra args
    assert ff_on_below(5,5,threshold=4,some_random_kwarg="test") == 0
    assert ff_on_below(3,5,threshold=4,some_random_kwarg="test") == 1


def test__get_ff_available():
    
    assert callable(FF_AVAILABLE["on"])
    assert FF_AVAILABLE["on"](1) == 1
    assert FF_AVAILABLE["on"](0) == 0
    
    assert callable(FF_AVAILABLE["off"])
    assert FF_AVAILABLE["off"](1) == 0
    assert FF_AVAILABLE["off"](0) == 1

    assert callable(FF_AVAILABLE["neutral"])
    assert FF_AVAILABLE["neutral"](1) == 1
    assert FF_AVAILABLE["neutral"](0) == 1

    assert callable(FF_AVAILABLE["on_above"])
    assert FF_AVAILABLE["on_above"](0.55,threshold=0.5) == 1
    assert FF_AVAILABLE["on_above"](0.45,threshold=0.5) == 0

    assert callable(FF_AVAILABLE["on_below"])
    assert FF_AVAILABLE["on_below"](0.55,threshold=0.5) == 0
    assert FF_AVAILABLE["on_below"](0.45,threshold=0.5) == 1
