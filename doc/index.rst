Welcome to scikit-guess
=======================

.. toctree::
   :maxdepth: 2
   :caption: Contents:


* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


Introduction
============

This scikit contains methods for computing fast, non-iterative estimates of
fitting parameters for common functions. The estimates may be used as-is on
their own, or refined through non-linear optimization algorithms. The name of
the scikit comes from the fact that estimates are a good initial guess for the
optimal fitting parameters.

The seed for this toolkit is Jean Jacquelin's paper "REGRESSIONS et EQUATIONS
INTEGRALES"\ [1]_. It demonstrates the techniques for deriving simple linear
least squares formulas by cleverly integrating the model functions. The
resulting estimates are not perfectly optimal, since the formulas are derived for
continuous functions but approximated by discrete samples. However, they are
very robust and fast. The key is to produce a result that is adequate for most
purposes, and can be used as a starting point for non-linear optimization
algorithms.


Installation
============

Installation through pypi is preferred, e.g. with ::

    pip install scikit-guess

Installation from source is also an option. The latest version can be obtained
from GitHub via ::

    git clone github.com:madphysicis/scikit-guess

You can then install in a couple of different ways. The standard is ::

    python setup.py install

To install for development purposes, or just to have the latest bleeding-edge
version always running ::

    python setup.py develop


Dependencies
============

This project strives to be a true scikit, and limit it's dependencies to just
numpy and scipy. Pandas support will be added at some point, but may not result
in a dependency, certainly not a mandatory one.

At this stage, the code is written in pure python, with all the extensions being
implemented through the dependencies. That may change at some point in the
future.

Python 2.7 and 3.4+ are supported.

Building the documentation requires sphinx.


Criteria
========

The algorithms presented in this scikit have two main criteria:

1. Fast
2. Non-iterative

Algorithms that provide a fast, non-iterative fit to functions that would
otherwise require non-linear fitting are welcome here. One required test is that
the function works faster than `scipy.optimize.curve_fit`.

While algorithms should yield good results, they do not need to be completely
optimal. The main criterion for accuracy is that the result can be used as a
good initial guess for a non-linear optimization method. As a consequece, some
of the algorithms provided by this scikit are good candidates for the `guess`
methods of corresponding models in `lmfit`.

All the functions that are *currently* supported are one dimensional, but that
is not a requirement.


Contributing
============

If you have an idea or would like to correct something, please submit an
`issue <https://github.com/madphysicist/scikit-guess/issues>`_, or a
`pull request <https://github.com/madphysicist/scikit-guess/pulls>`_.

If you have questions or are not sure if something is a bug, feel free to submit
an issue, or ask a question tagged ``scikit-guess`` on
`Stack Overflow <https://stackoverflow.com/>`_.

This scikit is still in early stages, so please provide all the criticism and
advice you can. Any support at all is welcome. In particular, the following
areas would be appreciated:

- pandas support
- Proper testing setup (e.g. TravisCI, Appveyor)
- Pointing out anything that is missing
- Adding new algorithms


Project Structure
=================

Each fitting algorithm resides in its own module. All the functions get imported
into the base `skg` namespace. Each module should contain a function called
`model` that applies the fitting parameters to a given set of x-values, either
raveled or along a particular axis (assuming the function is 1D). Multiple
algorithms that fit to the same model can live in the same module.


Testing
=======

Testing is semi-automated (still WIP). All the modules containing a fitting
function and a model will be tested against randomly generated inputs, and
checked for speed and quality. The quality of each algorithm will be assessed
based on these tests. Quality has three categories: speed, accuracy and
usefulness.

- Speed is a benchmark against `scipy.optimize.curve_fit`. An algorithm that is
  slower than a non-linear optimizer starting with default parameters is not
  very useful.
- Accuracy is checked by making sure that the fit is within reasonable bounds of
  the values computed by `scipy.optimize.curve_fit`. Reasonableness is a
  function of the analytically derived partial derivatives of the model with
  respect to the parameters.
- Usefulness is a measure of how many iterations `scipy.optimize.curve_fit`
  saves by using the algorithm as an initial guess. If the combined runtime of
  the algorithm and `curve_fit` is less than the runtime of `curve_fit` with
  default parameters, that's a win.


.. rubric:: Footnotes

.. [1] Jacquelin, Jean. "REGRESSIONS Et EQUATIONS INTEGRALES", pp. 15–18.,
   Available online https://www.scribd.com/doc/14674814/Regressions-et-equations-integrales
