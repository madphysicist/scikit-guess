"""
Shared utility functions used by the fitting routines.
"""

from numpy import (
    argsort, array, float_, inexact, issubdtype, __version__ as __np_version__
)
from numpy.core.multiarray import normalize_axis_index
from numpy.lib import NumpyVersion


__all__ = [
    'moveaxis', 'preprocess', 'preprocess_pair',
]


if NumpyVersion(__np_version__) >= '1.11.0':
    from numpy import moveaxis
else:
    from numpy import rollaxis
    def moveaxis(a, start, end):
        return rollaxis(a, start, normalize_axis_index(end, a.ndim + 1))


def preprocess(x, copy=False, float=False):
    """
    Ensure that `x` is a properly formatted numpy array.

    Proper formatting means at least one dimension, and may include
    optional copying and coersion into a floating point datatype.

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

    return array(x, copy=copy, subok=not copy, ndmin=1, dtype=dtype)


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
    """
    x = preprocess(x, copy=xcopy, float=True).ravel()
    y = preprocess(y, copy=ycopy, float=True).ravel()
    if x.shape != y.shape:
        raise ValueError('x and y must be the same shape')
    if not sorted:
        # Is there a better way to do this in scipy?
        ind = argsort(x)
        x = x[ind]
        y = y[ind]
    return x, y

