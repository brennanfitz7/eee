"""
Validate plot objects.
"""

from eee._private.interface import MockContextManager

import matplotlib
from matplotlib import pyplot as plt

def check_plot(fig=None,
               ax=None,
               context_manager=None):
    
    # Deal with fig, ax inputs, and widget inputs
    if fig is None and ax is None:
        fig, ax = plt.subplots(1,figsize=(6,6))
    
    if fig is not None and ax is None:
        err = "\n`ax` must be specified if `fig` is specified\n\n"
        raise ValueError(err)
    
    if not issubclass(type(ax),matplotlib.axes.Axes):
        err = "\n`ax` must be a matplotlib Axes object\n\n"
        raise ValueError(err)
    
    if fig is not None:
        if not issubclass(type(fig),matplotlib.figure.Figure):
            err = "\n`fig` must be a matplotlib Figure object\n\n"
            raise ValueError(err)

    if context_manager is None:
        context_manager = MockContextManager()

    if not hasattr(context_manager,"__enter__") or \
       not hasattr(context_manager,"__exit__"):
        err = "\n`context_manager` should be an object with both __enter__ and\n"
        err += "__exit__ methods\n\n"
        raise ValueError(err)

    return fig, ax, context_manager
