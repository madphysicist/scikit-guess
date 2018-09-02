"""
The namespace of scikit-guess is called `skg`. It exports all of the parameter
estimation functions provided by the scikit.

The functions themselves live in individual sub-modules.

.. autosummary::
   :toctree: generated/

   exp_fit
   pow_fit
"""

from .exp import exp_fit
from .pow import pow_fit


__all__ = ['exp_fit', 'pow_fit']

