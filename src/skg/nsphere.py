"""
N-sphere fit with center and radius.

This module is a little different from the others because it fits an
n-dimensional surface and because it does not have a model function
because of the non-functional nature of n-spheres.
"""

from __future__ import division, absolute_import

from numpy import asfarray, empty, sqrt, square
from numpy import __version__ as __np_version__
from numpy.lib import NumpyVersion
from scipy.linalg import lstsq


__all__ = ['nsphere_fit']


if NumpyVersion(__np_version__) >= '1.11.0':
    from numpy import moveaxis
else:
    from numpy import rollaxis
    def moveaxis(a, start, end):
        if end == -1:
            end = a.ndim
        elif end < 0:
            end += 1
        return rollaxis(a, start, end)


def nsphere_fit(x, axis=-1):
    r"""
    Fit an n-sphere to ND data.

    The center and radius of the n-sphere are optimized using the Coope
    method. The sphere is described by

    .. math::

       \left \lVert \vec{x} - \vec{c} \right \rVert_2 = r

    Parameters
    ----------
    x : array-like
        The n-vectors describing the data. Usually this will be a nxm
        array containing m n-dimensional data points.
    axis : int
        The axis that determines the number of dimensions of the
        n-sphere. All other axes are effectively raveled to obtain an
        ``(m, n)`` array.

    Return
    ------
    r : scalar
        The optimal radius of the best-fit n-sphere for `x`.
    c : array
        An array of size `x.shape[axis]` with the optimized center of
        the best-fit n-sphere.

    References
    ----------
    Coope, I.D. "CIRCLE FITTING BY LINEAR AND NONLINEAR LEAST SQUARES",
    Available online: https://core.ac.uk/download/pdf/35472611.pdf
    """
    x = asfarray(x)
    n = x.shape[axis]
    if axis not in (-1, x.ndim - 1):
        x = moveaxis(x, axis, -1)
    x = x.reshape(-1, n)
    m = x.shape[0]

    B = empty((m, n + 1), dtype=x.dtype)
    B[:, :-1] = x
    B[:, -1] = 1

    d = square(x).sum(axis=-1)

    y, *_ = lstsq(B, d, overwrite_a=True, overwrite_b=True)

    c = 0.5 * y[:-1]
    r = sqrt(y[-1] + square(c).sum())

    return r, c

