"""
Tests for the :func:`skg.gauss_fit` function.
"""

import numpy as np

from skg.gauss import gauss_fit


def test_paper():
    """
    Verifies the results of the example in :ref:`reei` in Section
    :ref:`reei1-sec3`.

    x- and y-values are copied from :ref:`reei-tab1`. Expected results
    are copied from :ref:`reei-fig1`.
    """
    x = [-0.992, -0.935, -0.836, -0.404, -0.326,
         -0.042,  0.068,  0.302,  0.439,  0.58]
    y = [0.238, 0.262, 0.38, 1.041, 0.922,
         0.755, 0.589, 0.34, 0.193, 0.083]

    mu_1, sigma_1 = -0.289356, 0.383915
    mu, sigma = gauss_fit(x, y)

    assert np.isclose(mu, mu_1, atol=5e-7, rtol=0.0)
    assert np.isclose(sigma, sigma_1, atol=5e-7, rtol=0.0)
