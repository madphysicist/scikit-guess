r"""
Sinusoidal fit based on the clustering algorithms in
:py:mod:`skg.cluster1d`.

The only parameter that is estimated with high accuracy is frequency.
All the other paramters can be obtained from other non-linear estimation
techniques.

.. math::

   f(x) = a sin(\omega t - \phi) + b

An alternative to this method is Jean Jacquelin's method.

.. todo::

   Add link to JJ's method above.

.. todo::

   Add proper handling of colinear inputs (and other singular matrix cases).

.. todo::

   Add proper handling for < 1 period.

.. todo::

   Add tests.

.. todo::

   Add nan_policy argument.

.. todo::

   Add axis parameter. Figure out how to do it properly.

.. todo::

   Add PEP8 check to formal tests.

.. todo::

   Allow broadcasting of x and y, not necessarily identical size
"""

from numpy import (
    argpartition, concatenate, full_like, searchsorted, median, std, sqrt,
    sin, pi
)

from .cluster1d import cluster1d
from .util import preprocess_pair, roots, counts, sums, medians, ends


__all__ = ['sin_fit']


def sin_fit(t, y, sorted=True, _debug=False):
    r"""
    Compute estimated initial parameters for noisy sinusoidal data.

    The result will have a very accurate approximation of the frequency,
    which is the crucial parameter for non-linear estimator initial
    values.

    Parameters
    ----------
    t : array-like
        The x-values of the data points. The fit will be performed on a
        raveled version of this array.
    y : array-like
        The y-values of the data points corresponding to `x`. Must be
        the same size as `x`. The fit will be performed on a raveled
        version of this array.
    sorted : bool
        Set to True if `x` is already monotonically increasing or
        decreasing. If False, `x` will be sorted into increasing order,
        and `y` will be sorted along with it.

    Return
    ------
    a : float
        Amplitude.
    omega : float
        Angular frequency: :math:`\omega = \frac{2\pi}{\lambda}`.
    phi : float
        Angular phase-shift.
    b : float
        Additive bias.

    References
    ----------
    - Currently none. This function, and the underlying clustering algorithm
      are entirely the work of the author. A peer reviewed paper is currently
      in the works.
    """
    # t: Time of raw readings
    # y: Sinusoidal readings
    t, y = preprocess_pair(t, y, sorted)
    # z: The interpolated crossings (in units of `t`)
    # bias: mean value of `y`
    z, bias = roots(t, y, return_bias=True)
    # breaks: indices into `z`, starting with zero
    breaks = cluster1d(z, sensitivity=2)
    # xings: medians of each group of `z` delimited by `breaks`, units of `t`
    xings = medians(z, breaks)

    if _debug:
        from matplotlib import pyplot as plt

        plt.figure()
        p1 = plt.subplot(3, 1, 1)
        p1.vlines(z, 0, 1)
        p1.vlines(0.5 * (z[breaks[1:]] + z[breaks[:-1] + 1]), 0, 1,
                   color='r', linestyle=':')
        p1.scatter(xings, full_like(xings, 0.5), c='m', zorder=10)
        p1.set_title('Initial grouping')

    # Post process crossings
    # Speed tested
    # indices: locations where `xings` fall into the original `t`
    indices = concatenate(([0], searchsorted(t, xings, side='right')))

    # sign: True for each region between a pair of `indices` that has at least
    #   half its points above the bias, False for the other regions. Same size
    #   as `indices` because of implicit last element.
    sign = (sums(y > bias, indices) > 0.5 * counts(indices, t.size))
    # mask: selection of regions that have opposite signs. Crossings that are
    #   between regions with the same sign are automatically deemed fake.
    mask = (sign[:-1] != sign[1:])

    if _debug:
        from matplotlib import pyplot as plt

        p2 = plt.subplot(3, 1, 2, sharex=p1)
        for sig, start, end in zip(sign, indices, ends(indices, t.size)):
            p2.plot(t[start:end], y[start:end], 'b' if sig else 'r')
        p2.plot(xings, full_like(xings, bias), 'mx', zorder=10)
        p2.plot(xings[~mask], full_like(xings[~mask], bias), 'mo', zorder=10)
        p2.set_title('Checking for sign flip')

    # Remove invalid `xings`
    xings = xings[mask]

    # Extract wavelengths
    wl_start = xings[:-2]
    wl_end = xings[2:]
    wavelengths = wl_end - wl_start
    wavelength = median(wavelengths)

    # Location of median wavelength    
    k = (wavelengths.size - 1) // 2
    wavelength_index = argpartition(wavelengths, k)[k]
    offset = (wavelength_index % 2 == 1) ^ sign[0]

    if _debug:
        p3 = plt.subplot(3, 1, 3, sharex=p1)
        for i, (start, end) in enumerate(zip(wl_start, wl_end)):
            color = 'g' if i == wavelength_index else 'r'
            p3.plot([start, end], [i % 3, i % 3], '-o', c=color)
        p3.set_title('Alternating wavelengths')

    omega = 2.0 * pi / wavelength
    amplitude = std(y) * sqrt(2.0) # Sine has crest factor ~= sqrt(2)
    full_phase = xings[wavelength_index + offset]
    phase = omega * (full_phase % wavelength)

    if _debug:
        lims = p2.get_xlim()
        p2.arrow(0, bias + 0.1 * amplitude, full_phase, 0, color='orange')
        p2.plot([full_phase] * 2, [bias + 0.12 * amplitude, 0], color='orange')
        p2.set_xlim(lims)

    return amplitude, omega, phase, bias


def model(t, a, omega, phi, b):
    r"""
    Compute :math:`a sin(\omega t - \phi) + b`.

    Parameters
    ----------
    t : array-like
        The value of the model will be the same shape as the input.
    a : float
        The amplitude of the sine wave at
        :math:`t = \frac{\left( \frac{\pi}{2} + \phi \right)}{\omega}`.
    omega : float
        The angular frequency of the sinusoid. This is :math:`2 \pi f`,
        where :math:`f` is the frequency.
    phi : float
        The angular phase shift of the sinusoid. This is the phase shift
        in units of periods. The phase shift can be expressed as
        :math:`t_0 = \frac{\phi}{\omega}`.
    b : float
        The bias of the sinusoid.

    Return
    ------
    y : array-like
        An array of the same shape as `t`, containing the model
        computed for the given parameters.
    """
    return a * sin(omega * t - phi) + b


sin_fit.model = model
