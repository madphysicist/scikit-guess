r"""
Weibull cumulative distribution fit with three parameters, of the form
:math:`F(t) = 1 - e^{-\left( \frac{t - \mu}{\beta} \right)^{\alpha}}`.

The data :math:`F_k` must be positive. Unlike the other estimations,
which are sorted according to their independent variable, this one
expects the dependent variable to be sorted in ascending order.


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

from numpy import array, expm1, log

from ._util import preprocess
from .exp import exp_fit


__all__ = ['weibull_cdf_fit']


def weibull_cdf_fit(x, y, sorted=True):
    r"""
    Weibull CDF fit of the form
    :math:`F(t) = 1 - e^{-\left( \frac{t - \mu}{\beta} \right)^{\alpha}}`.

    This implementation is based on a transformation that makes the
    approximate solution to integral equation :eq:`exp-eq`, presented in
    :ref:`ref-reei` applicable.

    Parameters
    ----------
    x : array-like
        The x-values of the data points, :math:`t` in the equation. The
        fit will be performed on a raveled version of this array.
    y : array-like
        The y-values of the data points corresponding to `x`,
        :math:`F(t)` in the equation. Must be the same size as `x`. The
        fit will be performed on a raveled version of this array.
    sorted : bool
        Set to True if `y` is already monotonically increasing or
        decreasing. If False, `y` will be sorted into increasing order,
        and `x` will be sorted along with it.

    Return
    ------
    alpha, beta, mu : ~numpy.ndarray
        A 3-element array of optimized fitting parameters. The first
        parameter is the shape, the second the scale, and the third the
        mean of the PDF.

    Notes
    -----
    The x- and y- data is conventionally transposed from what the final
    optimization expects. `y` is expected to be ascending or descending,
    not `x`.

    References
    ----------
    .. [1] Jacquelin, Jean. "\ :ref:`ref-reei`\ ",
       :ref:`pp. 19-20. <reei2-sec3>`,
       https://www.scribd.com/doc/14674814/Regressions-et-equations-integrales
    """
    F, y = preprocess(y, x, sorted)

    x = log(-log(1.0 - F))
    a, b, c = exp_fit(x, y, sorted=True)

    out = array([1.0 / c, b, a])
    return out


def model(x, alpha, beta, mu):
    r"""
    Compute :math:`y = 1 - e^{-\left( \frac{x - \mu}{\beta} \right)^{\alpha}}`.

    Parameters
    ----------
    x : array-like
        The value of the model will be the same shape as the input.
    alpha : float
        The shape of the distribution.
    beta : float
        The scale of the distribution.
    mu : float
        The mean of the distribution.

    Return
    ------
    y : array-like
        An array of the same shape as `x`, containing the model
        computed for the given parameters.
    """
    # TODO: Check if in-place operations are really faster here
    return -expm1(-((x - mu) / beta)**alpha)


weibull_cdf_fit.model = model
