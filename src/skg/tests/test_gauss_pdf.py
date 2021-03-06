"""
Tests for the :func:`skg.gauss_pdf_fit` function.
"""

import numpy as np

from skg.gauss_pdf import gauss_pdf_fit, model


def test_paper(plots):
    """
    Verifies the results of the example in :ref:`reei` in Section
    :ref:`reei1-sec3`.

    x/y-values and expected results are copied from
    :ref:`reei-gauss-pdf-data`. Results are displayed in
    :ref:`reei-gauss-pdf-plot`.
    """
    x = [-0.992, -0.935, -0.836, -0.404, -0.326,
         -0.042,  0.068,  0.302,  0.439,  0.58]
    y = [0.238, 0.262, 0.38, 1.041, 0.922,
         0.755, 0.589, 0.34, 0.193, 0.083]

    mu_1, sigma_1 = -0.289356, 0.383915
    mu, sigma = gauss_pdf_fit(x, y)

    # Parameters in the paper are rounded to 6 decimal places
    assert np.isclose(mu, mu_1, atol=5e-7, rtol=0.0)
    assert np.isclose(sigma, sigma_1, atol=5e-7, rtol=0.0)

    if plots:
        from .util import plotting_context, save

        x_2 = np.linspace(-1, 1, 100)
        with plotting_context() as fig:
            ax = fig.subplots(1)
            ax.scatter(x, y, color='k', marker='+')
            ax.plot(x_2, model(x_2, mu_1, sigma_1), c='k', ls='-')
            ax.plot(x_2, model(x_2, mu, sigma), c='r', ls=':')
            save(fig, __name__, 'paper')
