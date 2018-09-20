"""
Global configuration and fixture setup for pytest.
"""

from errno import EEXIST
from os import makedirs
from os.path import join
from warnings import warn

import numpy as np
from pytest import fixture

from .noise import white_normal, white_triangular, white_uniform


def pytest_addoption(parser):
    """
    Add some options to the default command line.

    The following options are added:

    `--plots`
        Draw plots of x-values, y-values and fit comparisons. This
        option checks if matplotlib is installed, and issues a warning
        if not.

    """
    parser.addoption("--plots", action="store_true", default=False,
                     help="Generate graphical plots of input data")


@fixture(scope='session')
def plots(request):
    """
    Enables debugging for the fixtures/tests that care about it.

    This fixture will only be set to `True` if the `--plots`
    command-line option is set.
    """
    if request.config.getoption('--plots'):
        try:
            makedirs(join(request.config.rootdir, '.skg_test'))
        except OSError as e:
            if e.errno != EEXIST:
                raise
        try:
            
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot
        except ImportError:
            warn('Unable to import matplotlib: '
                 'plotting will not be enabled.')
            return False

        return True

    return False


@fixture(scope='session')
def quality_metric():
    """
    Sets up the quality metric dictionary before running all tests, and
    publishes to an HTML file when complete.
    """
    quality_metric = {}
    yield quality_metric
    # TODO: Publish as HTML page


@fixture(scope='module', params=[
    0x0000, #0x1111, 0x1234, 0xBEEF, 0xCAFE, 0xDEAD, 0xFFFF
])
def seed(request):
    """
    Sets the seed for in numpy.random.

    Return the seed value, so it can be used by plots as part of the
    label.
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
