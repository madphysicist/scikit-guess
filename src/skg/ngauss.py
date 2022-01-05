r"""
N-dimensional, unnormalized Gaussian fit.

This method is adapded from the references to apply to any number of
dimensions, and arbitrary orientation of the error ellipsoid.

.. math::

   S(\vec{x}) = a e^{-\frac{1}{2} \left(\vec{x} - \vec{\mu}\right)^T \Sigma^{-1} \left(\vec{x} - \vec{\mu}\right)}

Here, :math:`\Sigma` is the covariance matrix that describes the
ellipsoid of the Gaussian. If this were a density function, the
amplitude ``a`` would be fixed at
:math:`\frac{1}{\sqrt{\left(2 \pi\right)^k det \Sigma}}`.

This function is specially suited for work with images. A wrapper to
facilitate image inputs is therefore provided:
:py:func:`ngauss_from_image`.
"""

from numpy import (
    add, einsum, empty, empty_like, exp, float_, indices, log, logical_not,
    reciprocal, square, std, triu_indices, zeros
)
from scipy.linalg import inv, lstsq

from .util import moveaxis, preprocess, preprocess_npair


__all__ = ['anthony_weights', 'ngauss_fit', 'ngauss_from_image']


def anthony_weights(noise=None):
    r"""
    Returns a function that computes weights for :py:func:`ngauss_fit`
    based on the `noise` standard deviation.

    Weights are computed (with appropriate clipping) according to

    .. math::

       w_i = \frac{1}{log(\frac{y_i + N}{y_i - N})}

    Signals smaller than `noise` are discarded by setting the weight
    to zero.

    Parameters
    ----------
    noise : float, optional
        An estimate of the standard deviation of signal noise. If not
        specified, use the standard deviation of `y`.

    Return
    ------
    weights : callable
        A function that accapts `x` and `y` and returns an array of
        weights the same size as `y`. `x` is completely ignored. This
        function is suitable as a `weights` parameter to
        :py:func:`ngauss_fit`.

    References
    ----------
    - [Anthony-Granick]_
    """
    def weights(x, y):
        # TODO: Time this.
        n = std(y) if noise is None else noise

        mask = y > n
        plus = y[mask]
        minus = plus.copy()
        plus += n
        minus -= n
        plus /= minus
        
        w = empty_like(y)
        log(plus, out=w, where=mask)
        reciprocal(w, out=w, where=mask)

        logical_not(mask, out=mask)
        w[mask] = 0
        return w


def ngauss_fit(x, y, axis=-1, weights=None, scaling=False):
    r"""
    An N-dimensional Gaussian fit.

    .. math::

        S(\vec{x}) = a e^{-\frac{1}{2} \left(\vec{x} - \vec{\mu}\right)^T \Sigma^{-1} \left(\vec{x} - \vec{\mu}\right)}

    This implementation is based on an extension of
    [AnthonyGranick]_ and the first of two passes described in
    [Wan-Wang-Wei-Li-Zhang]_.

    A weighing function must be applied to the data to avoid having the
    low-SNR data dominate the fit. The default is to weight the
    measurements by their intensity, as per _[Wan-Wang-Wei-Li-Zhang].
    However, other schemes are possible, such as the one proposed by
    _[Anthony-Granick]. The latter scheme can be used by passing in a
    callable returned by :py:func:`anthony_weights` as the `weights`
    parameter.

    Parameters
    ----------
    x : array-like
        The x-values of the data points. The fit will be performed on a
        version of the array raveled across all dimensions except
        `axis`.
    y : array-like
        The y-values of the data points corresponding to `x`. Must be
        the same size as `x` except at `axis`. The fit will be
        performed on a raveled version of this array.
    axis : int
        The axis containing the vectors of x. The dimension of the
        Gaussian is ``x.shape[axis]``. The default is ``-1``.
    weights : array-like or callable, optional
        Either an array with the same number of elements as `y` (it will
        be raveled), or a callable that accepts reshaped versions of `x`
        and `y` and returns an array of weights.
    scaling : bool, optional
        If `True`, scale and offset the data to a bounding box of -1 to
        +1 in each axis during computations for numerical stability.
        Default is `False`.

    Return
    ------
    a : float
        The amplitude of the Gaussian.
    mu : ~numpy.ndarray
        The mean of the Gaussian, as an N-element array.
    sigma : ~numpy.ndarray
        The covariance of the Gaussian, as an NxN positive definite
        matrix.

    See Also
    --------
    ngauss_from_image : A wrapper for fitting over an entire array.
    anthony_weights : A function that returns an alternative, noise-specific
        weighting scheme.

    Notes
    -----
    Negative and zero weights are discarded from the computation
    without ever being inserted into the solution matrix.

    References
    ----------
    - [Anthony-Granick]_
    - [Wan-Wang-Wei-Li-Zhang]_
    - :ref:`ngauss-suplement`, :ref:`ngauss-supplement-ndim`
    """
    x, y = preprocess_npair(x, y, axis, xcopy=False, ycopy=False)

    mask = y > 0
    x = x[mask, :]
    y = y[mask]

    if weights is None:
        weights = y
    elif callable(weights):
        weights = weights(x, y)

    mask = weights > 0
    x = x[mask, :]
    y = y[mask]
    weights = weights[mask]

    m, n = x.shape
    i = n * (n + 1) // 2  # Size of upper tri. of cov matrix

    if scaling:
        xmin = x.min(axis=0)
        xmax = x.max(axis=0)
        scale = 0.5 * (xmax - xmin)
        offset = 0.5 * (xmax + xmin)
        x -= offset
        x /= scale

    M = empty((m, i + n + 1), dtype=x.dtype)
    ind = triu_indices(n)
    M[:, :i] = x[:, ind[0]] * x[:, ind[1]]
    M[:, i:-1] = x
    M[:, -1] = 1
    M *= weights[:, None]

    p = log(y)
    p *= weights

    param, *_ = lstsq(M, p, overwrite_a=True, overwrite_b=True)

    sigma = zeros((n, n), dtype=param.dtype)
    sigma[ind] = -param[:i]
    add(sigma, sigma.T, out=sigma)
    sigma = inv(sigma, overwrite_a=True)

    mu = sigma @ param[i:-1]

    amp = exp(param[-1] + 0.5 * einsum('i,ij,j->', mu, sigma, mu))

    if scaling:
        mu *= scale
        mu += offset
        sigma *= scale
        sigma *= scale[:, None]

    return amp, mu, sigma


def ngauss_from_image(img, weights=None, scaling=True):
    """
    Compute a Gaussian fit to an entire image.

    Parameters
    ----------
    img : array-like
        The image to process. Usually a segment of a 2D image. The data
        is expected to have been background subtracted and thresholded
        so that any low-SNR pixels are set to zero.
    weights : array-like or callable, optional
        A weighing function must be applied to the data to avoid having
        the low-SNR data dominate the fit. The default is to weight
        the measurements by their intensity, as per
        [Wan-Wang-Wei-Li-Zhang]_. However, other schemes are possible,
        such as the one proposed by [Anthony-Granick]_. `weights` can
        be passed in as an array with the same number of elements as
        `y` (it will be raveled), or a callable that accepts reshaped
        versions of `x` and `y` and returns an array of weights.
    scaling : bool, optional
        If `True`, scale and offset the data to a bounding box of -1 to
        +1 in each axis during computations for numerical stability.
        Default is `True`.

    Return
    ------
    a : float
        The amplitude of the Gaussian.
    mu : ~numpy.ndarray
        The mean of the Gaussian, as an N-element array.
    sigma : ~numpy.ndarray
        The covariance of the Gaussian, as an NxN positive definite
        matrix.
    """
    index = moveaxis(indices(img.shape, sparse=False, dtype=float_), 0, -1)
    mask = img > 0
    img = img[mask]
    index = index[mask, :]
    return ngauss_fit(index, img, axis=-1, weights=weights, scaling=scaling)


def model(x, a, mu, sigma, axis=-1):
    r"""
    Compute :math:`y = a e^{-\frac{1}{2}\left(\frac{x - \mu}{\sigma}\right)^2}`.

    The number of dimensions, N, is determined from ``x.shape[axis]``.
    `mu` must be a vector of length N, and `sigma` must be an NxN
    positive-definite matrix.

    Parameters
    ----------
    x : array-like
        The value of the model will be the same shape as the input,
        with `axis` reduced.
    a : float
        The amplitude at :math:`\vec{x} = \vec{\mu}`.
    mu : array-like
        The location of the peak. Must be a an array of length
        ``x.shape[axis]``. May be a scalar if the location is the same
        value in all dimensions. This feature should only be used for
        a peak at zero.
    sigma : float
        The covariance matrix of the Gaussian. Must be a
        positive-definite matrix of shape
        ``(x.shape[axis], x.shape[axis])``.
    axis : int
        The axis corresponding to the dimension of ``x`` that contains
        the point vectors.

    Return
    ------
    y : array-like
        An array of the same shape as `x` with `axis` reduced,
        containing the model computed for the given parameters.

    Example
    -------
    To generate a 100px x 100px 2D image with a spot in the middle:

        >>> import numpy as np
        >>> from skg import ngauss_fit
        >>> x = np.indices((100, 100), dtype=float)
        >>> cov = np.array([[100., 40.], [40., 64.]])
        >>> img = ngauss_fit.model(x, 255, (50, 50), cov, axis=0)
    """
    x = preprocess(x, float=True, copy=True, axis=axis)
    x -= mu
    arg = einsum('...i,ij,...j', x, inv(sigma), x)
    square(arg, out=arg)
    arg *= -0.5
    return a * exp(arg)


ngauss_fit.model = model

