from future import division, absolute_import

from numpy import log, power

from .exp import exp_fit


__all__ = ['pow_fit']


def pow_fit(x, y, sorted=True):
    """
    Fit a power curve to raveled 1D data.

    The fitting parameters are comptued for:

    .. math::

       y = A + Bx^C

    Parameters
    ----------
    x : array-like
        The x-values of the data points. The fit will be performed on a
        raveled version of this array. All elements must be positive.
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
        the third is the power.

    Notes
    -----
    ``pow_fit(x, y, sorted)`` is equivalent to
    ``exp_fit(log(x), y, sorted)`` since

    .. math::

       A + Be^{Cx} = A + B(e^x)^C

    The algorithm and the theory behind it is presented in the paper
    below.

    References
    ----------
    Jacquelin, Jean. "REGRESSIONS Et EQUATIONS INTEGRALES", pp. 15–18.,
    Available online: https://www.scribd.com/doc/14674814/Regressions-et-equations-integrales
    """
    return exp_fit(log(x), y, sorted)


def model(x, a, b, c):
    """
    Compute

    .. math::

       y = A + Bx^C

    Parameters
    ----------
    x : array-like
        The value of the model will be the same shape as the input.
    a : float
        The additive bias.
    b : float
        The multiplicative bias.
    c : float
        The power.

    Return
    ------
    y : array-like
        An array of the same shape as ``x``, containing the model
        computed for the given parameters.
    """
    return a + b * power(x, c)

