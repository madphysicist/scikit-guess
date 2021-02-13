r"""
Unnormalized Gaussian bell curve fit.

The amplitude of this function is one of the fitting parameters, unlike
for the two-parameter PDF version.

.. math::

   f(x) = a e^{-\frac{1}{2} \left(\frac{x - \mu}{\sigma}\right)^2}

The third fitting parameter, :math:`a`, is the amplitude of the Gaussian
at :math:`x = \mu`. This is equivalent, up to a scaling factor, to
normalizing the area under the curve, as the PDF version does.

The conversion between amplitude :math:`a` and normalization :math:`A`
is given in :ref:`reei-supplement-gauss3` as

.. math::

   a = \frac{A}{\sigma \sqrt{2 \pi}}

For for the normalized (two parameter) Gaussian probability density
function, see :mod:`~skg.gauss_pdf`. For the CDF, see
:mod:`~skg.gauss_cdf`.

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

from numpy import array, cumsum, diff, empty, exp, sqrt
from scipy.linalg import lstsq

from .util import preprocess_pair


__all__ = ['gauss_fit']


def gauss_fit(x, y, sorted=True):
    r"""
    Gaussian bell curve fit of the form
    :math:`a e^{-\frac{1}{2}\left(\frac{x - \mu}{\sigma}\right)^2}`.

    This implementation is based on an extentsion the approximate
    solution to integral equation :eq:`gauss-pdf-eq`, presented in
    :ref:`ref-reei` and extended in :ref:`reei-supplement-extended`.

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
    a, mu, sigma : ~numpy.ndarray
        A three-element array containing the estimated amplitude, mean
        and standard deviation, in that order.

    References
    ----------
    - [Jacquelin]_ "\ :ref:`ref-reei`\ ", :ref:`pp. 6-8. <reei1-sec3>`
    - :ref:`reei-supplement-extended`, :ref:`reei-supplement-gauss3`
    """
    x, y = preprocess_pair(x, y, sorted)

    d = 0.5 * diff(x)
    xy = x * y

    # Did a timeit. This is the fastest way I could find to fill the matrix
    M = empty(xy.shape + (2,), dtype=xy.dtype)
    M[0, :] = 0
    cumsum((y[1:] + y[:-1]) * d, out=M[1:, 0])
    cumsum((xy[1:] + xy[:-1]) * d, out=M[1:, 1])

    Y = y - y[0]

    (A, B), *_ = lstsq(M, Y, overwrite_a=True, overwrite_b=True)

    mu, sigma = -A / B, sqrt(-1.0 / B)

    # Timeit shows that this is faster than a2 = model(x, 1.0, mu, sigma)
    m = exp(-0.5 * ((x - mu) / sigma)**2)
    amp = y.dot(m) / m.dot(m)

    out = array([amp, mu, sigma])

    return out


def model(x, a, mu, sigma):
    r"""
    Compute :math:`y = a e^{-\frac{1}{2}\left(\frac{x - \mu}{\sigma}\right)^2}`.

    Parameters
    ----------
    x : array-like
        The value of the model will be the same shape as the input.
    a : float
        The amplitude at :math:`x = \mu`.
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
    return a * exp(-0.5 * ((x - mu) / sigma)**2)


gauss_fit.model = model

