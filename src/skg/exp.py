"""
Exponential fit with additive bias.

.. todo::

   Add proper handling of colinear inputs (and other singular matrix cases).

.. todo::

   Add tests.

.. todo::

   Add nan_policy argument.

.. todo::

   Simplify matrix accumulation and solving to use least squares to setup
   the matrix for you.

.. todo::

   Add PEP8 check to formal tests.

.. todo::

   Add axis parameter.
"""

from __future__ import absolute_import, division

from numpy import cumsum, diff, empty, empty_like, exp, square
from numpy.linalg import solve

from ._util import preprocess


__all__ = ['exp_fit']


def exp_fit(x, y, sorted=True):
    """
    Exponential fit of the form :math:`A + Be^{Cx}`.

    Parameters
    ----------
    x : array-like
        The x-values of the data points. The fit will be performed on a
        raveled version of this array.
    y : array-like
        The y-values of the data points corresponding to `x`. Must be
        the same size as `x`. The fit will be performed on a raveled
        version of this array.
    sorted : bool
        Set to True if `x` is already monotonically increasing or
        decreasing. If False, `x` will be sorted into increasing order,
        and `y` will be sorted along with it.

    Return
    ------
    a, b, c : ~numpy.ndarray
        A 3-element array of optimized fitting parameters. The first
        element is the additive bias, the second the multiplicative, and
        the third the exponential.

    References
    ----------
    .. [1] Jacquelin, Jean. "REGRESSIONS Et EQUATIONS INTEGRALES", pp. 15-18.,
       https://www.scribd.com/doc/14674814/Regressions-et-equations-integrales
    """
    x, y = preprocess(x, y, sorted)

    s = empty_like(y)
    s[0] = 0
    s[1:] = cumsum(0.5 * (y[1:] + y[:-1]) * diff(x))

    xn = x - x[0]
    yn = y - y[0]

    # Recast this as a simple least-squares projection
    sx2 = square(xn).sum()
    sxs = (xn * s).sum()
    sys = (yn * s).sum()
    ss2 = square(s).sum()
    sxy = (xn * yn).sum()

    out = empty(3, dtype=float)

    _, out[2] = solve([[sx2, sxs], [sxs, ss2]], [sxy, sys])

    ex = exp(out[2] * x)

    # Recast this as a simple least-squares projection
    se1 = ex.sum()
    se2 = square(ex).sum()
    sy0 = y.sum()
    sye = (y * ex).sum()

    out[0], out[1] = solve([[x.size, se1], [se1, se2]], [sy0, sye])

    return out


def model(x, a, b, c):
    """
    Compute

    .. math::

       y = A + Be^{Cx}

    Parameters
    ----------
    x : array-like
        The value of the model will be the same shape as the input.
    a : float
        The additive bias.
    b : float
        The multiplicative bias.
    c : float
        The exponent.

    Return
    ------
    y : array-like
        An array of the same shape as ``x``, containing the model
        computed for the given parameters.
    """
    return a + b * exp(c * x)


exp_fit.model = model

