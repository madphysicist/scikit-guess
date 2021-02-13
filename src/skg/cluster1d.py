"""
A one-dimensional clustering algorithm that is useful for many
applications. It can be used to decode UART packets and seed iterative
minimizers for sinusoidal data.

The algorithm is very simple, relying on one sort, two diffs, and either
a max or partition of the data. The sort makes the asymptotic complexity
increase as O(n log(n)), but with a remarkably small load factor.

This module is not exactly a fitting rougine, although it is used by
:py:func:`~skg.sin_fit`, so it is not exported through the main
:py:mod:`skg` namespace.
"""

from numpy import (
    argmax, argpartition, argsort, array, concatenate, diff, log1p
)


__all__ = ['cluster1d', 'dscale_log']


def cluster1d(x, absolute=True, dscale=None, sensitivity=0):
    """
    Perform a 1-D clustering on the data.

    Data is assumed to be sorted, either increasing or decreasing
    monotonically. All data will be treated as a raveled 1-D array.

    Parameters
    ----------
    x : array-like
        The data points. Treated as a one-dimensional raveled array.
        Must be sorted.
    absolute : bool, optional
        Whether or not to introduce an absolute scaling term to the
        edge length computation. The term will determine whether cluster
        boundaries are subject to minute fluctuations in the data. It is
        generally a good idea to have this flag turned on (it is by
        default).
    dscale : callable, optional
        A function that accepts an array as input, to be applied
        *in-place* to the edges lengths before partitioning. This
        callable will adjust the location of the partition, and
        therefore what is considered to be a cluster boundary. Useful
        choices are `None` (no-op), `dscale_log`.
    sensitivity : int, optional
        The number non-maximal second-order edges to consider. With
        ``sensitivity=0`` or ``sensitivity=1``, only the largest gaps
        are considered. For ``sensitivity=2``, the second largest
        category of gaps will be included.

    Return
    ------
    q : numpy.ndarray
        An array of indices of the start of each cluster. The first
        element is always zero. The indices are always monotonically
        increasing.

    Notes
    -----
    `q` is intended to be used as the `index` argument to
    :py:meth:`numpy.ufunc.reduceat`. It can also be used as the input to
    :py:func:`numpy.split` by stripping off the leading element:
    ``q[1:]``. Use the functions in :py:mod:`skg.util` to slice and dice
    the indices further.

    Although funcionally identical, it is likely that
    ``sensitivity=0`` is slightly more efficient than ``sensitivity=1``.

    References
    ----------
    - Currently none. This function, is entirely the work of the author.
      A peer reviewed paper is currently in the works.
    """
    x = array(x, copy=False, subok=True).ravel()

    dx = diff(x)

    if dscale:
        dscale(dx)

    ind = argsort(dx)
    s = dx[ind]

    if absolute:
        # Speed tested
        ds = s.copy()
        ds[1:] -= s[:-1]
    else:
        ds = diff(s)

    m = argpartition(ds, -sensitivity)[-sensitivity:].min() \
        if sensitivity else argmax(ds)

    if not absolute:
        m += 1

    ind[m:] += 1
    ind[m:].sort()
    if m > 0:
        ind[m - 1] = 0
        return ind[m - 1:]
    return concatenate(([0], ind))


def dscale_log(x):
    """
    Perform an in-place log-scaling of the input, adjusted so that zeros
    are handled correctly.

    In practice this function computes

    .. math::

       log(x + 1)

    There is no return value, as the function operates in-place.

    Parameters
    ----------
    x : array-like
        The array to scale logarithmically.
    """
    log1p(x, out=x)
