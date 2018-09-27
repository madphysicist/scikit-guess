"""
Shared utility functions used by the fitting routines.

.. note::

   I suspect that some of these are already implemented in numpy in one
   form or another. If so, I would be more than happy to use those
   versions instead.
"""

from numpy import argsort, asfarray, swapaxes
from numpy.core.multiarray import normalize_axis_index


def preprocess(x, y, axis=None, sorted=True):
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
    axis : int or None
        The axis along which arrays will be reduced. If `None`, `x`
        and `y` are raveled.
    sorted : bool
        Set to True if `x` is already monotonically increasing or
        decreasing. If False, `x` will be sorted into increasing order,
        and `y` will be sorted along with it.

    Return
    ------
    x, y : ~numpy.ndarray
        Normalized versions of the inputs.
    axis : int
        Always an integer. If the input was `None`, returns zero.
    """
    x = asfarray(x)
    y = asfarray(y)
    if axis is None:
        x = x.ravel()
        y = y.ravel()
        axis = 0
    if x.shape != y.shape:
        raise ValueError('x and y must be the same shape')
    if not sorted:
        # Is there a better way to do this in scipy?
        ind = argsort(x, axis=axis)
        x = x[ind]
        y = y[ind]
    return x, y, axis


def replace_axis(arr, axis, size):
    """
    Returns a shape tuple of the same size as `arr.shape`, but with the
    element at `axis` replaced by `size`.

    Parameters
    ----------
    arr : ~numpy.ndarray
        The array to use for a shape reference.
    axis : int
        The axis to replace. Must be an integer, but may be negative.
    size : int
        The new size along `axis`.

    Return
    ------
    shape : tuple
        The udpated shape.
    """
    axis = normalize_axis_index(axis, arr.ndim)
    shape = arr.shape
    return shape[:axis] + (size,) + shape[axis + 1:]


def axis_slice(arr, axis, index):
    """
    Returns an index into `arr` with only `axis` fixed to `index`.

    The leading axes are effectively indexed with ``:``, and the
    trailing ones with ``...``.

    The axis of interest can be indexed by a scalar, a `slice` object
    and even `None` or :obj:`np.newaxis`. The latter case can be used to
    insert a new dimensions into an array.

    Parameters
    ----------
    arr : ~numpy.ndarray
        The array to modify.
    axis : int
        The axis to hold fixed.
    index : int or slice or None
        The index to hold `axis` at. Whatever this object is, it will be
        placed directly into the index tuple.

    Return
    ------
    index : tuple
        A slice index such that the element at `axis` is `index` and all
        the other axes are effectively indexed with ``:``.
    """
    axis = normalize_axis_index(axis, arr.ndim)
    ind = (slice(None),) * axis + (index, Ellipsis)
    return ind


def transpose_last(arr):
    """
    Return a view of `arr` with the last two axes swapped.

    This is suitable for use with ``@`` and :func:`matmul`.

    Parameter
    ---------
    arr : ~numpy.ndarray
        The array to view, interpreted as a stack of matrices.

    Return
    ------
    trans : ~numpy.ndarray
        A view of `arr` with the last two axes swapped. 1D arrays are
        passed through.
    """
    if arr.ndim >= 2:
        arr = swapaxes(arr, -1, -2)
    return arr
