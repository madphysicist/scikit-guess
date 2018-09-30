"""
Gaussian probability density fit.

The function in this module is normalized to one, as a PDF should be:

.. math::

   f(x) = \frac{1}{\sigma \sqrt{2 \pi}}
       e^{-\frac{1}{2} \left(\frac{x - \mu}{\sigma}\right)^2} \quad
   \int_{-\infty}^{+\infty} f(t)dt = 1

For for the unnormalized Gaussian (with an additional amplitude
parameter), see :mod:`gauss`. For the CDF, see :mod:`gauss_cdf`.

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

from __future__ import absolute_import, division

from numpy import array, cumsum, diff, empty, exp, pi, sqrt
from scipy.linalg import lstsq

from ._util import preprocess


__all__ = ['gauss_pdf_fit']


def gauss_pdf_fit(x, y, sorted=True):
    r"""
    Gaussian PDF fit of the form

    .. math::

       \frac{1}{\sigma \sqrt{2 \pi}}
           e^{-\frac{1}{2}\left(\frac{x - \mu}{\sigma}\right)^2}

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
    .. [1] Jacquelin, Jean. "\ :ref:`ref-reei`\ ",
       :ref:`pp. 6-8. <reei1-sec3>`,
       https://www.scribd.com/doc/14674814/Regressions-et-equations-integrales
    """
    x, y = preprocess(x, y, sorted)

    d = 0.5 * diff(x)
    xy = x * y

    a = empty(xy.shape + (2,), dtype=xy.dtype)
    a[0, :] = 0
    a[1:, 0] = cumsum((y[1:] + y[:-1]) * d)
    a[1:, 1] = cumsum((xy[1:] + xy[:-1]) * d)

    b = y - y[0]

    (A, B), *_ = lstsq(a, b, overwrite_a=True, overwrite_b=True)
    out = array([-A / B, sqrt(-1.0 / B)])

    return out


def model(x, mu, sigma):
    """
    Compute

    .. math::

       y = \frac{1}{\sigma \sqrt{2 \pi}}
           e^{-\frac{1}{2}\left(\frac{x - \mu}{\sigma}\right)^2}

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
        An array of the same shape as ``x``, containing the model
        computed for the given parameters.
    """
    return exp(-0.5 * ((x - mu) / sigma)**2) / (sigma * sqrt(2 * pi))


gauss_pdf_fit.model = model

