"""
Noise generators for the tests.

All the generators should have a signature that accepts amplitude and
optionally the size of the output array.

This module should not contain any attributes that are not modules,
dunder attributes, or noise generators.
"""

import numpy as np


__all__ = ['white_normal', 'white_triangular', 'white_uniform']


def white_normal(amplitude, size=None):
    """
    Normally distributed white noise.

    Amplitude is the 1-sigma value rather than a hard limit.
    """
    return np.random.normal(loc=0.0, scale=amplitude, size=size)


def white_triangular(amplitude, size=None):
    """
    White noise with a symmetric triangular distribution.
    """
    return np.random.triangular(
        left=-amplitude, mode=0.0, right=amplitude, size=size
    )


def white_uniform(amplitude, size=None):
    """
    White noise with a uniform distribution between the amplitude
    bounds.
    """
    return np.random.uniform(low=-amplitude, high=amplitude, size=size)
