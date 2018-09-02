============
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
resulting estimates are not perfectly optimal, since the formulas are derived
for continuous functions but approximated by discrete samples. However, they
are very robust and fast. The key is to produce a result that is adequate for
most purposes, and can be used as a starting point for non-linear optimization
algorithms.


--------
Criteria
--------

The algorithms presented in this scikit have two main criteria:

1. Fast
2. Non-iterative

Algorithms that provide a fast, non-iterative fit to functions that would
otherwise require non-linear fitting are welcome here. One required test is
that the function works faster than :func:`scipy.optimize.curve_fit`.

While algorithms should yield good results, they do not need to be completely
optimal. The main criterion for accuracy is that the result can be used as a
good initial guess for a non-linear optimization method. As a consequece, some
of the algorithms provided by this scikit are good candidates for the
:meth:`~lmfit.model.Model.guess` methods of corresponding models in
:doc:`lmfit <lmfit:index>`.

All the functions that are *currently* supported are one dimensional, but that
is not a requirement.


.. rubric:: Footnotes

.. [1] Jacquelin, Jean. "REGRESSIONS Et EQUATIONS INTEGRALES", pp. 15–18.,
   Available online https://www.scribd.com/doc/14674814/Regressions-et-equations-integrales