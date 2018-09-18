"""
Global fixture setup for pytest.
"""

from errno import EEXIST
from os import makedirs
from warnings import warn

import numpy as np
from pytest import fixture

from .noise import white_normal, white_triangular, white_uniform


@fixture(scope='module', params=[
    0x0000, 0x1111, 0x1234, 0xBEEF, 0xCAFE, 0xDEAD, 0xFFFF
], autouse=True)
def seed(request):
    """
    Sets the seed for in numpy.random.

    Return the seed value.
    """
    seed = request.param
    np.random.seed(seed)
    return seed


@fixture(scope='module', params=[
    white_normal, white_uniform, white_triangular
])
def noise_distribution(request):
    """
    A sequence of noise distribution types.

    All distribututions are implemented as functions that accept an
    amplitude and optionally the output size. The distributions are
    centered and symmetrical about zero.

    1. White Gaussian (normal) noise: amplitude is the standard
       deviation.
    2. Uniform noise: amplitude is the upper and lower bounds of the
       distribution.
    3. Triangular noise: amplitude is the upper and lower bounds of the
       distribution.
    """
    return request.param


@fixture(scope='session')
def debug():
    """
    Enables debugging for the fixtures/tests that care about it.

    .. todo::

       Enable this fixture only if the --debug command line option is
       specified.
    """
    try:
        makedirs('.skg_test')
    except OSError as e:
        if e.errno != EEXIST:
            raise
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot
    except ImportError:
        warn('Unable to import matplotlib: '
             'graphical debugging will be disabled.')

    return True


@fixture(scope='session')
def quality_metric():
    """
    Sets up the quality metric dictionary before running all tests, and
    publishes to an HTML file when complete.
    """
    quality_metric = {}
    yield quality_metric
    # TODO: Publish as HTML page
