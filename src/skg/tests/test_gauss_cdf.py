"""
Tests for the :func:`skg.gauss_pdf_fit` function.
"""

import numpy as np

from skg.gauss_cdf import gauss_cdf_fit


def test_paper():
    """
    Verifies the results of the example in :ref:`reei` in Section
    :ref:`reei1-appendix2`.

    x/y-values and expected results are copied from
    :ref:`reei-gauss-cdf-data`. Results are displayed in
    :ref:`reei-gauss-cdf-plot`.
    """
    x = [-0.914, -0.556, -0.49, -0.195, 0.019,
         0.045,  0.587,  0.764,  0.81, 0.884]
    y = [0.001, 0.017, 0.021, 0.097, 0.258,
         0.258, 0.704, 0.911, 0.911, 0.979]

    mu_1, sigma_1 = 0.266843, 0.374462
    mu, sigma = gauss_cdf_fit(x, y)

    assert np.isclose(mu, mu_1, atol=5e-7, rtol=0.0)
    assert np.isclose(sigma, sigma_1, atol=5e-7, rtol=0.0)
