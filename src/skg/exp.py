"""
Exponential fit with and additive bias of the form :math:`A + Be^{Cx}`.

.. todo::

   Add proper handling of colinear inputs (and other singular matrix cases).

.. todo::

   Add tests.

.. todo::

   Add nan_policy argument.

.. todo::

   Add PEP8 check to formal tests.

.. todo::

   Add axis parameter.
"""

from __future__ import absolute_import, division

from numpy import array, cumsum, diff, empty, exp, subtract
from scipy.linalg import lstsq

from ._util import preprocess


__all__ = ['exp_fit']


def exp_fit(x, y, sorted=True):
    r"""
    Exponential fit of the form :math:`A + Be^{Cx}`.

    This implementation is based on the approximate solution to integral
    equation :eq:`exp-eq`, presented in :ref:`ref-reei`.

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
    .. [1] Jacquelin, Jean. "\ :ref:`ref-reei`\ ",
       :ref:`pp. 15-18. <reei2-sec2>`,
       https://www.scribd.com/doc/14674814/Regressions-et-equations-integrales
    """
    x, y = preprocess(x, y, sorted)

    M = empty(y.shape + (2,), dtype=y.dtype)
    subtract(x, x[0], out=M[:, 0])
    M[0, 1] = 0
    cumsum(0.5 * diff(x) * (y[1:] + y[:-1]), out=M[1:, 1])

    Y = y - y[0]

    (A, B), *_ = lstsq(M, Y, overwrite_a=True, overwrite_b=True)

    a, c = -A / B, B

    M[:, 0].fill(1.0)
    exp(c * x, out=M[:, 1])

    (a, b), *_ = lstsq(M, y, overwrite_a=True, overwrite_b=False)

    out = array([a, b, c])

    return out


def model(x, a, b, c):
    """
    Compute :math:`y = A + Be^{Cx}`.

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
        An array of the same shape as `x`, containing the model
        computed for the given parameters.
    """
    return a + b * exp(c * x)


exp_fit.model = model

