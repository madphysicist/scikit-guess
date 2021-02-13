"""
The namespace of scikit-guess is called :mod:`skg`. It exports all of
the parameter estimation functions provided by the scikit.

The functions themselves live in individual sub-modules. Each module
containing a fitting function also defines a `model` attribute, which
refers to the model function being fitted. Each fitting routine has a
`model` attribute pointing to the module-level `model`.

.. rubric:: Sub-modules

.. autosummary::
   :toctree: generated/

   exp
   gauss_cdf
   gauss_pdf
   gauss
   pow
   weibull_cdf
   nsphere

.. rubric:: Exported Functions

.. autosummary::
   :toctree: generated/

   exp_fit
   gauss_cdf_fit
   gauss_pdf_fit
   gauss_fit
   pow_fit
   weibull_cdf_fit
   nsphere_fit
   sin_fit
"""

from .exp import exp_fit as exp_fit
from .gauss_cdf import gauss_cdf_fit as gauss_cdf_fit
from .gauss_pdf import gauss_pdf_fit as gauss_pdf_fit
from .gauss import gauss_fit as gauss_fit
from .pow import pow_fit as pow_fit
from .weibull_cdf import weibull_cdf_fit as weibull_cdf_fit
from .nsphere import nsphere_fit as nsphere_fit

from .version import __version__ as __version__


__all__ = [
    'gauss_cdf_fit', 'gauss_pdf_fit', 'gauss_fit', 'exp_fit', 'pow_fit',
    'weibull_cdf_fit', 'nsphere_fit',
]


def test(*args, **kwargs):
    """
    Run the tests.

    Positional arguments will be inserted as command line arguments to
    the main test routine. Keyword arguments will be passed directly.
    """
    from pytest import main

    cmd = ['-p', 'skg.tests.options', '--pyargs', 'skg.tests']
    cmd.extend(args)
    return main(cmd, **kwargs)
