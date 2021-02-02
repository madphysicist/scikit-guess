"""
Tests for the :func:`skg.nsphere_fit` function.
"""

import numpy as np

from skg.nsphere import nsphere_fit


def test_paper():
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

    # 8-points, no scaling
    r, c = nsphere_fit(x[:8], axis=-1, scaling=False)
    assert np.isclose(r, r8, atol=5e-6, rtol=0.0)
    assert np.allclose(c, c8, atol=5e-6, rtol=0.0)

    # 9-points, no scaling
    r, c = nsphere_fit(x, axis=-1, scaling=False)
    assert np.isclose(r, r9, atol=5e-6, rtol=0.0)
    assert np.allclose(c, c9, atol=5e-6, rtol=0.0)

    # 8-points, with scaling
    r, c = nsphere_fit(x[:8], axis=-1, scaling=True)
    assert np.isclose(r, r8, atol=5e-6, rtol=0.0)
    assert np.allclose(c, c8, atol=5e-6, rtol=0.0)

    # 9-points, with scaling
    r, c = nsphere_fit(x, axis=-1, scaling=True)
    assert np.isclose(r, r9, atol=5e-6, rtol=0.0)
    assert np.allclose(c, c9, atol=5e-6, rtol=0.0)

