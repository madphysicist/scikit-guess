"""
Shared utility functions used by the fitting routines.
"""

from numpy import argsort, asfarray, swapaxes
from numpy.core.multiarray import normalize_axis_index


def preprocess(x, y, sorted=True):
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

    Return
    ------
    x, y : ~numpy.ndarray
        Normalized versions of the inputs.
    """
    x = asfarray(x).ravel()
    y = asfarray(y).ravel()
    if x.shape != y.shape:
        raise ValueError('x and y must be the same shape')
    if not sorted:
        # Is there a better way to do this in scipy?
        ind = argsort(x)
        x = x[ind]
        y = y[ind]
    return x, y

