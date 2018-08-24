"""
Exponential fit with additive bias.

.. todo::

   Add proper handling of colinear inputs (and other singular matrix cases).

.. todo::

   Add tests.

.. todo::

   Add nan_policy argument.
"""

from future import absolute_import, division

from numpy import (
    argsort, asfarray, cumsum, diff, empty, empty_like, exp, square,
)
from numpy.linalg import inv


__all__ = ['exp_fit']


def exp_fit(x, y, sorted=True):
    """
    Fit an exponential curve to raveled 1D data.

    This algorithm does not require any a-priori knowledge of the data,
    such as the intercept. The fitting parameters are comptued for:

    .. math::

       y = A + Be^{Cx}

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
    a, b, c : array
        A 3-element array of optimized fitting parameters. The first
        element is the additive bias, the second the multiplicative, and
        the third the exponential.

    References
    ----------
    Jacquelin, Jean. "REGRESSIONS Et EQUATIONS INTEGRALES", pp. 15–18.,
    Available online: https://www.scribd.com/doc/14674814/Regressions-et-equations-integrales
    """
    x = asfarray(x).ravel()
    y = asfarray(y).ravel()
    if x.size != y.size:
        raise ValueError('x and y must be the same size')
    if not sorted:
        # Is there a better way to do this in scipy?
        ind = argsort(x)
        x = x[ind]
        y = y[ind]

    s = empty_like(y)
    s[0] = 0
    s[1:] = cumsum(0.5 * (y[1:] + y[:-1]) * diff(x))

    xn = x - x[0]
    yn = y - y[0]

    sx2 = square(xn).sum()
    sxs = (xn * s).sum()
    sys = (yn * s).sum()
    ss2 = square(s).sum()
    sxy = (xn * yn).sum()

    out = empty(3, dtype=float)

    _, out[2] = inv([[sx2, sxs], [sxs, ss2]]).dot([[sxy], [sys]])

    ex = exp(out[2] * x)

    se1 = ex.sum()
    se2 = square(ex).sum()
    sy0 = y.sum()
    sye = (y * ex).sum()

    out[0], out[1] = inv([[x.size, se1], [se1, se2]]).dot([[sy0], [sye]])

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

