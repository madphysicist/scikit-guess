"""
Shared utility functions used by the fitting routines.
"""

from numpy import (
    add, argsort, array, concatenate, empty, flatnonzero, float_, inexact,
    issubdtype, mean, median, sqrt, __version__ as __np_version__
)
from numpy.core.multiarray import normalize_axis_index
from numpy.lib import NumpyVersion


__all__ = [
    'bias_mean', 'counts', 'ends', 'medians', 'moveaxis', 'preprocess',
    'preprocess_pair', 'preprocess_npair', 'roots', 'stds', 'sums', 'vars'
]


if NumpyVersion(__np_version__) >= '1.11.0':
    from numpy import moveaxis
else:
    from numpy import rollaxis
    def moveaxis(a, start, end):
        return rollaxis(a, start, normalize_axis_index(end, a.ndim + 1))


def preprocess(x, copy=False, float=False, axis=None):
    """
    Ensure that `x` is a properly formatted numpy array.

    Proper formatting means at least one dimension, and may include
    optional copying, reshaping and coersion into a floating point
    datatype.

    Parameters
    ----------
    x : array-like
        The array to process. If not already a numpy array, it will be
        converted to one.
    copy : bool, optional
        If True, a copy is made regardless of whether `x` is already a
        numpy array or not. The default is False.
    float : bool, optional
        If True, and `x` is not an inexact array already
        (:py:attr:`numpy.float16`, :py:attr:`numpy.float32`,
        :py:attr:`numpy.float64`, :py:attr:`numpy.float96`,
        :py:attr:`numpy.float128`, etc), coerce to be of type
        :py:attr:`numpy.float_`. Defaults to False.
    axis : int, optional
        If specified, the specified axis is moved to the end of the
        shape. Default is to return `x` with the original dimensions.

    Return
    ------
    x : ~numpy.ndarray
        Processed version of the input.
    """
    if float:
        dtype = x.dtype if hasattr(x, 'dtype') and \
                           issubdtype(x.dtype, inexact) else float_
    else:
        dtype=None

    x = array(x, copy=copy, subok=not copy, ndmin=1, dtype=dtype)
    if axis is not None and axis not in (-1, x.ndim - 1):
        # moveaxis always returns a new view, never the same object
        x = moveaxis(x, axis, -1)
    return x


def preprocess_pair(x, y, sorted=True, xcopy=False, ycopy=False):
    """
    Ensure that `x` and `y` are floating point arrays of the same size,
    ranked in increasing order by `x`.

    Parameters
    ----------
    x : array-like
        The x-values of the data points. The array will be converted to
        floating point, raveled and sorted, only as necessary.
    y : array-like
        The y-values of the data points corresponding to `x`. Must be
        the same size as `x`. Will be converted to floating point and
        raveled only as necessary. Will be sorted if `x` gets sorted.
    sorted : bool
        Set to True if `x` is already monotonically increasing or
        decreasing. If False, `x` will be sorted into increasing order,
        and `y` will be sorted along with it.
    xcopy : bool, optional
        Ensure that `x` gets copied even if it is already an array. The
        default is to leave arrays untouched as much as possible.
    ycopy : bool
        Ensure that `y` gets copied even if it is already an array. The
        default is to leave arrays untouched as much as possible.

    Return
    ------
    x, y : ~numpy.ndarray
        Processed versions of the inputs.

    See Also
    --------
    preprocess_npair : Similar function but for `x` containing vectors
        and `y` scalars.
    """
    x = preprocess(x, copy=xcopy, float=True)
    y = preprocess(y, copy=ycopy, float=True)
    if x.shape != y.shape:
        raise ValueError('x and y must be the same shape')
    x = x.ravel()
    y = y.ravel()
    if not sorted:
        # Is there a better way to do this in scipy?
        ind = argsort(x)
        x = x[ind]
        y = y[ind]
    return x, y


def preprocess_npair(x, y, axis=-1, xcopy=False, ycopy=False):
    """
    Ensure that `x` and `y` are floating point arrays of compatible
    size.

    `x` is an array containing vectors along dimension `axis`. `y`
    contains scalar elements. The shape of `y` must match that of `x`
    exactly except for `axis`.

    Parameters
    ----------
    x : array-like
        The vector x-values of the data points. The array will be
        converted to floating point, and raveled along all dimensions
        but `axis`, which will be the last dimension.
    y : array-like
        The y-values of the data points corresponding to `x`. Must have
        one fewer dimension than `x`, and its shape must match all
        elements of `x`'s shape except `axis`. Will be converted to
        floating point and raveled.
    xcopy : bool, optional
        Ensure that `x` gets copied even if it is already an array. The
        default is to leave arrays untouched as much as possible.
    ycopy : bool
        Ensure that `y` gets copied even if it is already an array. The
        default is to leave arrays untouched as much as possible.

    Return
    ------
    x, y : ~numpy.ndarray
        Processed versions of the inputs.

    See Also
    --------
    preprocess_pair : For cases when `x` and `y` both contain scalars,
        and are the exact same size.
    """
    if axis is None:
        raise ValueError('Axis must be an integer, not None')
    x = preprocess(x, copy=xcopy, float=True, axis=axis)
    y = preprocess(y, copy=ycopy, float=True)
    if x.shape[:-1] != y.shape:
        raise ValueError('x and y must be the same shape besides axis in x')
    x = x.reshape(-1, x.shape[-1])
    y = y.ravel()
    return x, y


def ends(p, size):
    """
    Compute the exclusive end indices of each cluster for a dataset of
    `size` elements.

    Parameters
    ----------
    p : array-like
        The start points of each cluster. Clusters continue until the
        start of the next cluster.
    size : int
        The size of the dataset (not of `p`). Necessary because `p`
        omits the end index.

    Return
    ------
    q : numpy.ndarray
        The exclusive indices of the end of each cluster.
        ``q.size == p.size``.
    """
    return concatenate((p.ravel()[1:], [size]))


def counts(p, size):
    """
    Compute the number of elements in each cluster in a dataset of `size`
    elements.

    Parameters
    ----------
    p : array-like
        The start points of each cluster. Clusters continue until the
        start of the next cluster.
    size : int
        The size of the dataset (not of `p`). Necessary because `p`
        omits the end index.

    Return
    ------
    c : numpy.ndarray
        The number of elements in each cluster. ``c.size == p.size``.
    """
    c = ends(p, size)
    c -= p.ravel()
    return c


def sums(x, p):
    """
    Compute the sums of each cluster in `x`.

    Parameters
    ----------
    x : array-like
        The array to compute sums over.
    p : array-like
        The start points of each cluster. Clusters continue until the
        start of the next cluster. The indices should start with zero
        and increase monotonically.

    Return
    ------
    s : numpy.ndarray
        The sum of each cluster identified by `p`.
    """
    return add.reduceat(preprocess(x).ravel(), p.ravel())


def means(x, p):
    """
    Compute the means of each cluster in `x`.

    Parameters
    ----------
    x : array-like
        The array to compute means over.
    p : array-like
        The start points of each cluster. Clusters continue until the
        start of the next cluster. The indices should start with zero
        and increase monotonically.

    Return
    ------
    m : numpy.ndarray
        The mean of each cluster identified by `p`.
    """
    x = preprocess(x)
    return sums(x, p) / counts(p, x.size)


def vars(x, p):
    """
    Compute the variances of each cluster in `x`.

    Parameters
    ----------
    x : array-like
        The array to compute variances over.
    p : array-like
        The start points of each cluster. Clusters continue until the
        start of the next cluster. The indices should start with zero
        and increase monotonically.

    Return
    ------
    v : numpy.ndarray
        The variance of each cluster identified by `p`.
    """
    x = preprocess(x)
    c = counts(p, x.size)
    return (sums(x**2, p) - sums(x, p)**2 / c) / c


def stds(x, p):
    """
    Compute the standard deviations of each cluster in `x`.

    Parameters
    ----------
    x : array-like
        The array to compute standard deviations over.
    p : array-like
        The start points of each cluster. Clusters continue until the
        start of the next cluster. The indices should start with zero
        and increase monotonically.

    Return
    ------
    s : numpy.ndarray
        The standard deviation of each cluster identified by `p`.
    """
    return sqrt(vars(x, p))


def medians(x, p):
    """
    Compute the median of each cluster in `x`.

    Parameters
    ----------
    x : array-like
        The array to compute medians over.
    p : array-like
        The start points of each cluster. Clusters continue until the
        start of the next cluster. The indices should start with zero
        and increase monotonically.

    Return
    ------
    m : numpy.ndarray
        The median of each cluster identified by `p`.
    """
    x = preprocess(x).ravel()
    m = empty(p.size)
    for i, (start, end) in enumerate(zip(p, ends(p, x.size))):
        m[i] = median(x[start:end])
    return m


def bias_mean(x, y):
    """
    Computes the bias of a dataset as the mean of the `y` values.

    This function is a wrapper to be used as `bias` for `roots`.

    Parameters
    ----------
    x : array-like
        The x-values, passed in to fulfill the correct interface and
        otherwise ignored.
    y : array-like
        The y-values of the data, treated as a 1D raveled array.

    Return
    ------
    The mean of `y`.
    """
    return mean(y)


def roots(x, y, sorted=True, bias=bias_mean, return_indices=False, return_bias=False):
    """
    Interpolate the roots of a 1-D dataset.

    Roots are interpolated using linear interpolation about an arbitrary
    bias, possibly generated from the data.

    Parameters
    ----------
    x : array-like
        The x-values of the data. Must be monotonically increasing or
        decreasing. Treated as a 1D raveled array.
    y : array-like
        The y-values of the data. Treated as a 1D raveled array. Must be
        the same lendth as `x`.
    sorted : bool
        Set to True if `x` is already monotonically increasing or
        decreasing. If False, `x` will be sorted into increasing order,
        and `y` will be sorted along with it.
    bias : scalar or array-like or callable, optional
        Either a fixed y-value that the data is offset by, or a callable
        that generates the value from `x` and `y`. The roots are
        computed for the y-values with the bias subtracted off. Defaults
        to the mean of the data.
    return_indices : bool, optional
        If `True`, return the indices of the interpolation points.
    return_bias : bool, optional
        If `True`, return the bias as an additional output parameter.
        This parameter is especially useful if `bias` is a callable.

    Return
    ------
    roots : numpy.ndarray
        The x-values of the interpolated intersections between the data
        and the bias.
    indices : numpy.ndarray
        If `return_indices` is set, this will be the indices of the
        interpolation points in `x`, as if returned by
        `numpy.searchsorted`. This will be the index of the point on the
        right of the interval containing the root.
    bias : numpy.ndarray or scalar
        If `return_bias` is set, this will be the actual scalar or array
        that is subtraced from the y-values to compute the roots.
    """
    if callable(bias):
        bias = bias(x, y)

    x, y = preprocess_pair(x, y, sorted, ycopy=True)
    # Speed tested
    y -= bias

    # Speed tested
    le = (y <= 0)
    gt = ~le

    mask = (gt[:-1] & le[1:]) | (le[:-1] & gt[1:])

    # Speed tested
    m0 = flatnonzero(mask)
    m1 = m0 + 1

    x0 = x[m0]
    x1 = x[m1]
    y0 = y[m0]
    y1 = y[m1]

    roots = (x0 * y1 - x1 * y0) / (y1 - y0)

    if return_indices:
        if return_bias:
            return roots, m1, bias
        return roots, m1
    if return_bias:
        return roots, bias
    return roots
