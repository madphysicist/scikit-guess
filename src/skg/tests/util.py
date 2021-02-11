"""
Collection of testing utilities.
"""

from contextlib import contextmanager


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
