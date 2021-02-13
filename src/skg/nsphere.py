"""
N-sphere fit with center and radius.

This module is a little different from the others because it fits an
n-dimensional surface and because it does not have a model function
because of the non-functional nature of n-spheres.
"""

from numpy import empty, sqrt, square
from numpy import __version__ as __np_version__
from numpy.lib import NumpyVersion
from scipy.linalg import lstsq

from .util import preprocess

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


def nsphere_fit(x, axis=-1, scaling=False):
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
    scaling : bool
        If `True`, scale and offset the data to a bounding box of -1 to
        +1 during computations for numerical stability. Default is
        `False`.

    Return
    ------
    r : scalar
        The optimal radius of the best-fit n-sphere for `x`.
    c : array
        An array of size `x.shape[axis]` with the optimized center of
        the best-fit n-sphere.

    References
    ----------
    - [Coope]_ "\ :ref:`ref-cfblanls`\ "
    """
    x = preprocess(x, float=True)
    n = x.shape[axis]
    if axis not in (-1, x.ndim - 1):
        x = moveaxis(x, axis, -1)
    x = x.reshape(-1, n)
    m = x.shape[0]

    B = empty((m, n + 1), dtype=x.dtype)
    X = B[:, :-1]
    X[:] = x
    B[:, -1] = 1

    if scaling:
        xmin = X.min()
        xmax = X.max()
        scale = 0.5 * (xmax - xmin)
        offset = 0.5 * (xmax + xmin)
        X -= offset
        X /= scale

    d = square(X).sum(axis=-1)

    y, *_ = lstsq(B, d, overwrite_a=True, overwrite_b=True)

    c = 0.5 * y[:-1]
    r = sqrt(y[-1] + square(c).sum())

    if scaling:
        r *= scale
        c *= scale
        c += offset

    return r, c

