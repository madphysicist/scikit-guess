"""
:mod:`skg.tests` is the `scikit-guess` test package. Each module in the main
:mod:`skg` package has a corresponding test module here.

Each function is tested for speed and quality using a modified version of
:func:`scipy.optimize.curve_fit` as the gold standard. The reason for choosing
:func:`~scipy.optimize.curve_fit` is that it is a very common and accessible
function to use. If this package provides benefits to users of
:func:`~scipy.optimize.curve_fit`, it will likely improve the experience for
other packages like `lmfit`_ as well.

Each test computes a quality metric that is stored in a global dictionary keyed
by test ID. The dictionary is set up by the
:func:`tests.conftest.quality_metric` fixture. Quality metrics are aggregated
and published in an HTML file upon completion.

The following tests are run for every fitting function in the scikit:

Speed
    A simple benchmark of the fit against the modified
    :func:`~scipy.optimize.curve_fit`.
Convergence
    Checks that :func:`~scipy.optimize.curve_fit` converges faster using
    the fit as an initial guess than it would with default parameters.
Total Speed
    Benchmarks the fit + :func:`~scipy.optimize.curve_fit` with the fit
    as an initial guess against just :func:`~scipy.optimize.curve_fit`
    with the default guess. Passing this test is optional.
RMS
    The RMS of the data about the function fitted with the fitting
    routine must be within an acceptable threshold of the RMS of the
    data about the reqult of :func:`~scipy.optimize.curve_fit`.
Accuracy
    The parameters computed by the fitting function must be within
    acceptable thresholds of the parameters computed by
    :func:`~scipy.optimize.curve_fit`.
Pathological Cases
    Each function has pathological cases that it can not handle
    properly. Tests for each function should be included on a
    case-by-case basis to ensure the contractal behavior in these cases.
    For example, exponential fits can not be made to colinear data,
    unless it is horizontal.
CurveFit
    An implicit test of :func:`~scipy.optimize.curve_fit` is done with
    each dataset, to ensure that the RMS of the data about the fit is
    lower than the RMS of the data about the model with initial fitting
    parameters (only for noisy datasets).
Input Parameters
    Each input parameter should be tested individually.

If the source material on wich a the algorithm is based provides a sample, a
"Paper Test" may be added to verify the results against what should be an
independent implementation.


.. include:: /link-defs.rst
"""
