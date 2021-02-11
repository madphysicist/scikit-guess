"""
Tests for the :func:`skg.nsphere_fit` function.
"""

import numpy as np

from skg.nsphere import nsphere_fit

from .util import plotting_context


def test_paper(plots):
    """
    Verifies the results of the example in the paper.

    x/y-values are found on p. 2 and expected results are on p. 4.
    Both the 8-point and 9-point cases are tested in this method.
    Tests are run with and without scaling.
    """

    # Last point is extra
    x = [[0.7, 4.0], [3.3, 4.7], [5.6, 4.0], [7.5, 1.3], [6.4, -1.1],
         [4.4, -3.0], [0.3, -2.5], [-1.1, 1.3], [3.0, 1.0]]
    r8, c8 = 4.10914, (3.06030, 0.74361)
    r9, c9 = 3.87132, (3.10253, 0.75467)

    x = np.array(x)

    # 8-points, no scaling
    r8ns, c8ns = nsphere_fit(x[:8], axis=-1, scaling=False)
    # Parameters in the paper are rounded to 5 decimal places
    assert np.isclose(r8ns, r8, atol=5e-6, rtol=0.0)
    assert np.allclose(c8ns, c8, atol=5e-6, rtol=0.0)

    # 9-points, no scaling
    r9ns, c9ns = nsphere_fit(x, axis=-1, scaling=False)
    # Parameters in the paper are rounded to 5 decimal places
    assert np.isclose(r9ns, r9, atol=5e-6, rtol=0.0)
    assert np.allclose(c9ns, c9, atol=5e-6, rtol=0.0)

    # 8-points, with scaling
    r8s, c8s = nsphere_fit(x[:8], axis=-1, scaling=True)
    # Parameters in the paper are rounded to 5 decimal places
    assert np.isclose(r8s, r8, atol=5e-6, rtol=0.0)
    assert np.allclose(c8s, c8, atol=5e-6, rtol=0.0)

    # 9-points, with scaling
    r9s, c9s = nsphere_fit(x, axis=-1, scaling=True)
    # Parameters in the paper are rounded to 5 decimal places
    assert np.isclose(r9s, r9, atol=5e-6, rtol=0.0)
    assert np.allclose(c9s, c9, atol=5e-6, rtol=0.0)

    if plots:
        from matplotlib.patches import Circle

        with plotting_context() as fig:
            ax = fig.subplots(2, 2)

            ax[0, 0].scatter(x[:8, 0], x[:8, 1], color='k', marker='x')
            ax[0, 0].add_patch(Circle(c8, r8, ec='k', fc='none', ls='-'))
            ax[0, 0].add_patch(Circle(c8ns, r8ns, ec='r', fc='none', ls=':'))
            ax[0, 0].set_title('8-Points, No Scaling')

            ax[0, 1].scatter(x[:, 0], x[:, 1], color='k', marker='x')
            ax[0, 1].add_patch(Circle(c9, r9, ec='k', fc='none', ls='-'))
            ax[0, 1].add_patch(Circle(c9ns, r9ns, ec='r', fc='none', ls=':'))
            ax[0, 1].set_title('9-Points, No Scaling')

            ax[1, 0].scatter(x[:8, 0], x[:8, 1], color='k', marker='x')
            ax[1, 0].add_patch(Circle(c8, r8, ec='k', fc='none', ls='-'))
            ax[1, 0].add_patch(Circle(c8s, r8s, ec='r', fc='none', ls=':'))
            ax[1, 0].set_title('8-Points, Scaling')

            ax[1, 1].scatter(x[:, 0], x[:, 1], color='k', marker='x')
            ax[1, 1].add_patch(Circle(c9, r9, ec='k', fc='none', ls='-'))
            ax[1, 1].add_patch(Circle(c9s, r9s, ec='r', fc='none', ls=':'))
            ax[1, 1].set_title('9-Points, Scaling')

            fig.savefig('.skg_test/{}-paper.debug.png'.format(__name__), dpi=300)

