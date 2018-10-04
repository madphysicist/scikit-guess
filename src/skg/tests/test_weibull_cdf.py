"""
Tests for the :func:`skg.weibull_cdf_fit` function.
"""

import numpy as np

from skg.weibull_cdf import weibull_cdf_fit


def test_paper():
    """
    Verifies the results of the example in :ref:`reei` in Section
    :ref:`reei2-sec3`.

    x/y-values and expected results are copied from
    :ref:`reei-weibull-cdf-data`. Results are displayed in
    :ref:`reei-weibull-cdf-plot`.
    """
    x = [1.202, 1.397, 1.537, 1.57, 1.768, 1.856, 1.87, 1.889, 1.918, 2.098,
         2.212, 2.349, 2.453, 2.557, 2.596, 2.602, 2.678, 2.706, 3.089, 3.441]
    y = [0.033, 0.082, 0.131, 0.181, 0.23, 0.279, 0.328, 0.377, 0.426, 0.475,
         0.524, 0.573, 0.622, 0.671, 0.72, 0.769, 0.818, 0.867, 0.916, 0.965]

    alpha_1, beta_1, mu_1 = 2.44301, 1.55262, 0.82099
    alpha, beta, mu = weibull_cdf_fit(x, y)

    assert np.isclose(alpha, alpha_1, atol=5e-6, rtol=0.0)
    assert np.isclose(beta, beta_1, atol=5e-6, rtol=0.0)
    assert np.isclose(mu, mu_1, atol=5e-6, rtol=0.0)
