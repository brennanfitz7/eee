
__authors__ = "Brennan Fitzgerald, Michael J. Harms"

from . import ddg
from . import core
from . import structure
from . import io
from . import calcs
from . import analysis
from . import utilities

def _check_for_notebook():
    """
    Check whether the code is being executed in a notebook or standard
    standard python interpreter.

    Return
    ------
        string for jupyter or IPython, None for something not recognized.
    """

    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            return "jupyter"   # Jupyter notebook or qtconsole
        elif shell == 'TerminalInteractiveShell':
            return "IPython"   # Terminal running IPython
        else:
            return None        # Not sure what interpreter

    # Probably standard Python interpreter
    except NameError:
        return None

_in_notebook = _check_for_notebook()

from .__version__ import __version__


