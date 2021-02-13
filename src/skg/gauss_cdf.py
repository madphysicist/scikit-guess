r"""
Gaussian cumulative distribution fit.

The function in this module is asymptotic to zero at negative infinity
and to one at positive infinity, as a CDF should be:

.. math::

   f(x) & = \frac{1}{\sqrt{2 \pi} \sigma} \int_{-\inf}^x exp\left(
            -\frac{1}{2}\left(\frac{t - \mu}{\sigma}\right)^2\right)dt \\
        & = \frac{1}{2} +
            \frac{1}{2} erf \left( \frac{x - \mu}{\sqrt{2} \sigma} \right)

A fit to the probability density function for this cumulative
distribution is provided in :mod:`~skg.gauss_cdf`. For for the
unnormalized Gaussian bell curve (with an additional amplitude
parameter), see :mod:`~skg.gauss`.

.. todo::

   Add proper handling of colinear inputs (and other singular matrix cases).

.. todo::

   Add tests.

.. todo::

   Add nan_policy argument.

.. todo::

   Add axis parameter. Figure out how to do it properly.

.. todo::

   Add PEP8 check to formal tests.

.. todo::

   Include amplitude in integrals.

.. todo::

   Allow broadcasting of x and y, not necessarily identical size
"""

from numpy import array, ones_like, sqrt, stack
from scipy.linalg import lstsq
from scipy.special import erf, erfinv

from .util import preprocess_pair


__all__ = ['gauss_cdf_fit']


def gauss_cdf_fit(x, y, sorted=True):
    r"""
    Gaussian CDF fit of the form
    :math:`\frac{1}{2}+\frac{1}{2}erf\left(\frac{x-\mu}{\sqrt{2}\sigma}\right)`.

    This implementation is based on the approximate solution to integral
    equation :eq:`gauss-cdf-fx2`, presented in :ref:`ref-reei`.

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
    mu, sigma : ~numpy.ndarray
        A two-element array containing the estimated mean and standard
        deviation, in that order.

    References
    ----------
    - [Jacquelin]_ "\ :ref:`ref-reei`\ ", :ref:`pp. 11-13. <reei1-appendix2>`
    """
    x, y = preprocess_pair(x, y, sorted)

    M = stack((x, ones_like(x)), axis=1)
    Y = erfinv(2 * y - 1)

    (A, B), *_ = lstsq(M, Y, overwrite_a=True, overwrite_b=True)
    out = array([-B / A, 1 / (sqrt(2.0) * A)])

    return out


def model(x, mu, sigma):
    r"""
    Compute :math:`y = \frac{1}{2} + \frac{1}{2} erf \left( \frac{x - \mu}{\sqrt{2} \sigma} \right)`.

    Parameters
    ----------
    x : array-like
        The value of the model will be the same shape as the input.
    mu : float
        The mean.
    sigma : float
        The standard deviation.

    Return
    ------
    y : array-like
        An array of the same shape as `x`, containing the model
        computed for the given parameters.
    """
    return 0.5 * (1 + erf((x - mu) / (sqrt(2) * sigma)))


gauss_cdf_fit.model = model

