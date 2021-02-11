r"""
Power fit with additive bias of the form :math:`A + Bx^C`.

As a general rule, ``pow_fit(x, y, ...)`` is equivalent to
``exp_fit(log(x), y, ...)`` since
:math:`A + Be^{Cx} = A + B \left( e^x \right)^C`.


.. todo::

   Add proper handling of colinear inputs (and other singular matrix cases).

.. todo::

   Add tests.

.. todo::

   Add `nan_policy` argument.
"""

from numpy import log, power

from .exp import exp_fit


__all__ = ['pow_fit']


def pow_fit(x, y, sorted=True):
    r"""
    Power fit of the form :math:`A + Bx^C`.

    This implementation is based on the approximate solution to integral
    equation :eq:`exp-eq`, presented in :ref:`ref-reei`. A power fit is
    regarded as an exponential fit with a logarithmically scaled x-axis
    in this algorightm.

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
    a, b, c : ~numpy.ndarray
        A three-element array containing the estimated additive and
        multiplicative biases and power, in that order.

    References
    ----------
    .. [1] Jacquelin, Jean. "\ :ref:`ref-reei`\ ",
       :ref:`pp. 15-18. <reei2-sec2>`,
       https://www.scribd.com/doc/14674814/Regressions-et-equations-integrales
    """
    return exp_fit(log(x), y, sorted)


def model(x, a, b, c):
    """
    Compute :math:`y = A + Bx^C`.

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
        An array of the same shape as `x`, containing the model
        computed for the given parameters.
    """
    return a + b * power(x, c)


pow_fit.model = model

