"""
The namespace of scikit-guess is called `skg`. It exports all of the
parameter estimation functions provided by the scikit.

The functions themselves live in individual sub-modules. Each module
containing a fitting function also defines a `model` attribute, which
refers to the model function being fitted. Each fitting routine has a
`model` attribute pointing to the module-level `model`.

.. autosummary::
   :toctree: generated/

   exp_fit
   pow_fit
"""

from .exp import exp_fit
from .pow import pow_fit

from .version import __version__


__all__ = ['exp_fit', 'pow_fit']

