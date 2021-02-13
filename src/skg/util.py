"""
Shared utility functions used by the fitting routines.
"""

from numpy import (
    argsort, array, float_, inexact, issubdtype
)


__all__ = [
    'preprocess', 'preprocess_pair',
]


def preprocess(x, copy=False, float=False):
    if float:
        dtype = x.dtype if hasattr(x, 'dtype') and \
                           issubdtype(x.dtype, inexact) else float_
    else:
        dtype=None

    return array(x, copy=copy, subok=not copy, ndmin=1, dtype=dtype)


def preprocess_pair(x, y, sorted=True, xcopy=False, ycopy=False):
    """
    Ensures that `x` and `y` are floating point arrays of the same size,
    ranked in increasing order by `x`.

    Parameters
    ----------
    x : array-like
        The x-values of the data points. The array will be converted to
        floating point, raveled and sorted, only as necessary.
    y : array-like
        The y-values of the data points corresponding to `x`. Must be
        the same size as `x`. Will be converted to floating point and
        raveled only as necessary. Will be sorted if `x` is sorted.
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
        Normalized versions of the inputs.
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

