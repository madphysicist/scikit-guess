"""
Collection of testing utilities for common actions.
"""

from contextlib import contextmanager
from os.path import join


FOLDER = '.skg_test'


@contextmanager
def plotting_context(*args, **kwargs):
    """
    Context manager that imports matplotlib as necessary, creates a
    figure, and closes it on exit.

    All parameters are passed directly to
    :py:func:`~matplotlib.pyplot.figure`.
    """
    from matplotlib import pyplot as plt
    figure = plt.figure(*args, **kwargs)
    yield figure
    plt.close(figure)


def save(fig, module, label, debug=False, dpi=300):
    """
    Saves a figure in the standard test output directory.

    The figure name is nomalized with the name of the calling module.
    """
    ext = 'debug.png' if debug else 'png'
    label = f'{module}-{label}.{ext}'

    fig.savefig(join(FOLDER, label), dpi=dpi)
