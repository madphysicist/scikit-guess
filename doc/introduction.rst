============
Introduction
============

This scikit contains methods for computing fast, non-iterative estimates of
fitting parameters for common functions. The estimates may be used as-is on
their own, or refined through non-linear optimization algorithms. The name of
the scikit comes from the fact that estimates are a good initial guess for the
optimal fitting parameters.

The seed for this toolkit is Jean Jacquelin's paper "Régressions et équations
intégrales"\ :ref:`[ref] <ref-reei>`. It demonstrates the techniques for
deriving linear least squares formulas by cleverly integrating the model
functions. The resulting estimates are not always optimal, since the
formulas are derived by discretizing continuous functions. However, they are
very robust and fast. The key is to produce a result that is adequate for most
purposes, and can be used as a starting point for non-linear optimization
algorithms.


--------
Criteria
--------

The algorithms presented in this scikit have two main criteria:

1. Fast
2. Non-iterative

Algorithms that provide a fast, non-iterative fit to functions that would
otherwise require non-linear fitting are welcome here. One required test is
that the function works faster than a modified version of
:func:`scipy.optimize.curve_fit`, with default initial parameters.

While algorithms should yield good results, they do not need to be completely
optimal. The main criterion for accuracy is that the result can be used as a
good initial guess for a non-linear optimization method. Some parameters are
more significant than others (think frequency of a sine wave). As a consequece,
some of the algorithms provided by this scikit are good candidates for the
:py:meth:`~lmfit.model.Model.guess` methods of corresponding models in
`lmfit`_.


---------------
Historical Note
---------------

I came across Jean Jacquelin's paper while researching fast exponential
fitting routines that could be run on enormous dataframes within a limited time
slot: i.e., stock applications. The subsequent attempt to submit
`PR #9158 <https://github.com/scipy/scipy/pull/9158>`_ to scipy_ was rejected.
The maintainers kindly indicated that the correct place for such a function
would be in an independent scikit rather the main body of scipy. At the same
time, Matt Newville, the author of lmfit_, noticed the PR, and we briefly
discussed using the functions of this scikit as the
:py:meth:`~lmfit.model.Model.guess` methods for the appropriate models. Ever
since then one of the goals of this package has been to get into a shape that
is useable by lmfit_.


.. include:: /link-defs.rst
