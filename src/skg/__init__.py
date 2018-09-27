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
   gauss_fit
   pow_fit
"""

from .gauss import gauss_fit
from .exp import exp_fit
from .pow import pow_fit

from .version import __version__


__all__ = ['gauss_fit', 'exp_fit', 'pow_fit']


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
