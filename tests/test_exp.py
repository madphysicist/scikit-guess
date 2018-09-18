"""
Tests for the :func:`skg.exp_fit` function.
"""

import sys

import numpy as np
from pytest import fixture

from skg.exp import exp_fit


partials = {
    'a': lambda x, *param: np.ones_like(x),
    'b': lambda x, *param: np.exp(param[2] * x),
    'c': lambda x, *param: x * np.exp(param[2] * x),
}


@fixture(scope='module', params=[10**x for x in range(1, 6)])
def n_points(request):
    """
    Number of points in the test data set.
    """
    return request.param


@fixture(scope='module', params=[0.0, 0.5, 1.0, 2.0])
def x_spread(request):
    """
    Ratio of :math:`\\mu` to :math:`\\sigma` of gamma distribution.

    The spread uniquely determines the shape since

    .. math::

       \\frac{\\mu}{\\sigma} = \\frac{\\theta k}{\\theta \\sqrt{k}} = \\frac{1}{\\sqrt{k}}
    """
    return request.param


@fixture(scope='module', params=[10**x for x in range(4)])
def x_range(request):
    """
    The range on the number line for x-value.

    This tests scaling. For exponentials, the range is -scale to +scale.
    """
    return -request.param, request.param


@fixture(scope='module')
def x_data(seed, n_points, x_spread, x_range, debug):
    """
    Generate x-data points either uniformly or with random spacing
    across a given range.

    If `x_spread` is zero, create a :func:`~np.linspace`. Otherwise,
    create a sequence of steps that progress with a
    :func:`~np.random.gamma` distribution.

    Generated sequences are logged to an image in the ``.skg_test``
    folder.
    """
    start, end = x_range
    if x_spread:
        # Gamma spacing
        d = end - start
        spread2 = x_spread**2
        space = np.random.gamma(1.0 / spread2, spread2 * d / n_points, size=n_points)
        space = np.cumsum(space) - space[0]
        x = start + space * d / space[-1]
    else:
        # Uniform spacing
        x = np.linspace(start, end, n_points)

    if debug and 'matplotlib.pyplot' in sys.modules:
        from matplotlib import pyplot as plt
        fig, ax = plt.subplots()
        ax.plot(range(x.size), x)
        label = 'Gamma with spread {}'.format(x_spread) if x_spread else 'Uniform'
        ax.set_title('X-VALUES\n{}, {} points\nFrom {} to {}'.format(label, n_points, start, end))
        fig.savefig('.skg_test/x_data_N{}_S{}_{}-{}.debug.png'.format(n_points, x_spread, start, end))
        plt.close(fig)

    return x


@fixture(scope='module')
def fitting_params(seed):
    # Generate params based on what?
    a = np.random.normal(scale=125.0) # Really doesn't matter
    b = np.random.normal(scale=25.0)  # Somewhat matters, but not much
    c = np.random.normal(scale=5.0)   # Matters the most
    return a, b, c


@fixture(scope='module')
def clean_data(x_data, fitting_params):
    return exp_fit.model(x_data, *fitting_params), fitting_params


@fixture(scope='module', params=[0.1, 0.25, 0.5, 1.0])
def noisy_data(request, seed, clean_data, noise_distribution):
    clean_data, fitting_params = clean_data
    amplitude = request.param * np.ptp(clean_data)
    noise = noise_distribution(amplitude, clean_data.shape)
    return clean_data + noise, fitting_params


@fixture(scope='module')
def horizontal_data(x_data):
    value = np.random.normal(scale=10.0)
    data = np.empty_like(x_data)
    data.fill(value)
    return data, (value, 0.0, 0.0)


@fixture(scope='module')
def colinear_data(x_data):
    m = np.random.normal(scale=0.1)
    b = np.random.normal(scale=10.0)
    return m * x_data + b, (np.nan, np.nan, np.nan)


def test_exp_clean(x_data, clean_data):
    y_data, fitting_param = clean_data
    print('OK')


def test_exp_noisy(x_data, noisy_data):
    y_data, fitting_param = noisy_data
    print('OK')
