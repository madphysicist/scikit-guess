"""
N-sphere fit with center and radius.

This module is a little different from the others because it fits an
n-dimensional surface and because it does not have a model function
because of the non-functional nature of n-spheres.
"""

from __future__ import division, absolute_import

from numpy import asfarray, sqrt
from numpy import __version__ as __np_version__
from numpy.lib import NumpyVersion
from scipy.linalg import lstsq, norm


__all__ = ['nsphere_fit']


if NumpyVersion(__np_version__) >= '1.11.0':
    from numpy import moveaxis
else:
    from numpy import rollaxis as moveaxis


def nsphere_fit(x, axis=0):
    r"""
    Fit an n-sphere to ND data.

    The center and radius of the n-sphere are optimized using the Coope
    method. The sphere is described by

    .. math::

       \|\vec{x} - \vec{c}\|_2 = r

    Parameters
    ----------
    x : array-like
        The n-vectors describing the data. Usually this will be a nxm
        array containing m n-dimensional data points.
    axis : int
        The axis that determines the number of dimensions of the
        n-sphere. All other axes are effectively raveled to obtain an
        nxm array.

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
    if axis:
        x = moveaxis(x, axis, 0)
    x = x.reshape(n, -1)
    m = x.shape[1]

    B = np.array((m, n + 1), dtype=x.dtype)
    B[:, :-1] = x.T
    B[:, -1] = 1

    D = np.square(x).sum(axis=0)

    y, *_ = lstsq(B, D, overwrite_a=True, overwrite_b=True)

    c = 0.5 * y[:-1]
    r = sqrt(y[-1] + norm(c, axis=0))

    return r, c

