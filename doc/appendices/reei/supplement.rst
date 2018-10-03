=======================
Supplementary Materials
=======================

.. contents::
   :local:

-------------------------
Additional Considerations
-------------------------

The purpose of this translation was originally to motivate the writing of the
scikit to which it is attached. As such, the calculations shown for linear
regressions in the paper can be greatly simplified using the power of existing
linear algebra solvers. The particular implementations discussed in this
section are targeted towards `numpy`_ and `scipy`_ in particular. However, the
formulation shown here will most likely be helpful for a number of other modern
linear least-squares packages. This is especially true since many software
packages, including `numpy`_ and `scipy`_, share `LAPACK`_ as a backend for
performing the actual computations.

With a linear algebra package such as `LAPACK`_, it is sufficient to compute
the coefficient matrix :math:`A` and ordinate vector :math:`b`, without
performing any further operations on them. A linear regression is generally
solved as something equivalent to

.. math::

   x = \left( M^T * M \right)^{-1} \left( M^T * p \right)

The formulae in the paper show how to compute the elements of :math:`M^T * M`
and :math:`M^T * p`, which is usually done more efficiently by exiting packages
when given the raw :math:`M` and :math:`p`. These variable names were chosen to
avoid conflict with the names :math:`A`, :math:`B`, :math:`a`, :math:`b`, etc,
which are used heavily by the paper and software package documentation for
different things.

The following sections show how to construct such simplified solutions to the
equations in the paper. Solutions are described briefly, and presented
concretely with Python code. The solutions here are for conceptual reference
only. They are not a complete or exact reflection of how things are done in the
scikit itself.

In the code snippets below, ``x`` and ``y`` are assumed to be one-dimensional
numpy arrays. Both arrays have an equal number of elements, ``n``. Numpy is
implicitly imported under the conventional name ``np``.

Gaussian PDF
============

Algorithm originally presented in :ref:`reei1-sec3` and summarized
:ref:`here <reei1-sec3-alg>`.

:math:`M` is a matrix with the cumulative sums :math:`S` and :math:`T` as
the columns. In numpy terms:

.. code-block:: python

   S = np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(x))
   S = np.insert(S, 0, 0)

   xy = x * y
   T = np.cumsum(0.5 * (xy[1:] + xy[:-1]) * np.diff(x))
   T = np.insert(T, 0, 0)

   M = np.stack((S, T), axis=1)

:math:`p` is the vector of measured :math:`y` values decreased by its first
element. In numpy terms:

.. code-block:: python

   p = y - y[0]

Gaussian CDF
============

Algorithm originally presented in :ref:`reei1-appendix2` and summarized
:ref:`here <reei1-appendix2-alg>`.

:math:`M` is a matrix with just the :math:`x` values and ones as the columns.
In numpy terms:

.. code-block:: python

   A = np.stack((x, np.ones_like(x)), axis=1)

:math:`p` is a more complicated function of :math:`y` in this case:

.. code-block:: python

   p = scipy.special.erfinv(2 * y - 1)

Exponential
===========

Algorithm originally presented in :ref:`reei2-sec2` and summarized
:ref:`here <reei2-sec2-alg>`.

In the first regression to obtain the parameter :math:`c`, we have for
:math:`M`:

.. code-block:: python

   S = np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(x))
   S = np.insert(S, 0, 0)

   M1 = np.stack((x - x[0], S), axis=1)

:math:`p` is the vector of measured :math:`y` values decreased by its first
element. In numpy terms:

.. code-block:: python

   p1 = y - y[0]

Once we find a fixed value for :math:`c_1`, we can fit the remaining parameters
with :math:`M` constructed from :math:`c_1` and :math:`x`. In numpy terms:

.. code-block:: python

   M2 = np.stack((np.ones_like(x), np.exp(c1 * x)), axis=1)

:math:`p` is just the raw :math:`y` values for this regression. In numpy terms:

.. code-block:: python

   p2 = y

.. include:: ../page_break.rst


.. _reei-supplement-extended:

---------------------
Extended Applications
---------------------

Some additional common functions and modifications that can be added to the
suite presented in the paper are presented here.


.. _reei-supplement-gauss3:

Three-Parameter Gaussian
========================

In the real world, data often comes out as an unnormalized Gaussian of three,
rather than two parameters, with the third parameter being the amplitude or
normalization. The normalization constant of the Gaussian PDF is chosen so that
the total area under the curve is unity:

.. math::

   \frac{1}{\sigma \sqrt{2 \pi}} \int_{-\infty}^{+\infty}
       e^{-\frac{1}{2} \left(\frac{t - \mu}{\sigma} \right)^2}dt = 1

The amplitude of the PDF at the peak :math:`x = \mu` is therefore
:math:`\frac{1}{\sigma \sqrt{2 \pi}}`. An unnormalized Gaussian can be
conveniently parametrized either in terms of its total area :math:`A`, or peak
amplitude :math:`a`, in addition to :math:`\mu` and :math:`\sigma`. We have the
relationship:

.. math::

   a = \frac{A}{\sigma \sqrt{2 \pi}}

Either is a valid choice for the third parameter. The scikit uses the amplitude
:math:`a` because it is easier to visualize, and makes the computation just a
little easier. We therefore wish to use the methodology presented in the
original paper to find the parameters :math:`a`, :math:`\mu` and :math:`\sigma`
which optimize the fit to some set of data points
:math:`(x_1, y_1), (x_2, y_2), ..., (x_k, y_k), ..., (x_n, y_n)` of our
function:

.. math::

   f(x) = a e^{-\frac{1}{2} \left(\frac{x - \mu}{\sigma} \right)^2}

The formulation of the integral equation :eq:`gauss-pdf-eq` is still
applicable because the amplitude parameter is not resolvable from those
equations. We are able to solve for the approximations :math:`\mu_1` and
:math:`\sigma_1` as before (but with a different :math:`f(x)`, which absorbs
the multiplicative constant).

Optimizing :math:`a` will require an additional step. Fortunately, with
:math:`\mu` and :math:`\sigma` fixed, the minimization of residuals in terms of
:math:`a` is already a linear problem:

.. math::

   \sum_{k=1}^n \varepsilon_k^2 = \sum_{k=1}^n \left(y_k - a E_k)\right)^2 \\
   \text{where} \quad E_k =
       e^{-\frac{1}{2} \left(\frac{x_k - \mu_1}{\sigma_1} \right)^2}

Setting the partial derivative of the sum of squared errors with respect to
:math:`a` to zero:

.. math::

   \sum_{k=1}^n \left( y_k - a E_k \right) E_k = 0

Since we are solving for only one unknown, a single equation is produced:

.. math::

   a = \frac{\sum_{k=1}^n y_k E_k}{\sum_{k=1}^n E_k^2}

.. include:: ../page_break.rst


------
Errata
------

The following list of errata was complied during the translation. All items
link back to the corresponding location in the translated paper as footnotes.

.. [errata-reei-1] The original list of function conventions read

   .. math::

      \Phi_k = \Phi(x_k); G_k = G(x_k); H_k = H(x_k); ...;
          \alpha_k = \alpha(x_k); \beta_k = \beta(x_k); ...

   The portion :math:`; \alpha_k = \alpha(x_k); \beta_k = \beta(x_k); ...` does
   not appear to be correct. The list is referring to functions of :math:`x_k`,
   which :math:`\alpha` and :math:`\beta` are definitely not. Subscripting them
   is probably an error.

.. [errata-reei-2] The original equations :eq:`gauss-pdf-int` and :eq:`gauss-pdf-eq`
   read:

   .. math::
      :label: gauss-pdf-int-old

      \int_{x_1}^x \left(t - \mu\right)f(t)dt =
          -\sqrt{\frac{\pi}{2}}\sigma\left(f(x) - f(x_1)\right)

   .. math::
      :label: gauss-pdf-eq-old

      \begin{cases}
          f(x) - f(x_1) = A \int_{x_1}^x f(t)dt + B \int_{x_1}^x t f(t)dt \\
          \text{with:} \quad A = \frac{\mu}{\sigma}\sqrt{\frac{2}{\pi}} \quad
          \text{and} \quad B = -\frac{1}{\sigma}\sqrt{\frac{2}{\pi}}
      \end{cases}

   Taking the integral :math:`\int_{x_1}^x \left(t - \mu\right) f(t)dt`, with
   :math:`f(x)` given in :eq:`gauss-pdf-fx`, shows that the factor of
   :math:`\sqrt{\frac{\pi}{2}}` should actually be an additional factor of
   :math:`\sigma`. This factor then propagates to :math:`A` and :math:`B` as
   shown in the translated paper.

.. [errata-reei-3] The original equation :eq:`gauss-pdf-lsq` had
   :math:`y - y_1` in the rightmost matrix:

   .. math::
      :label: gauss-pdf-lsq-old

      \begin{bmatrix}
          \sum \left(y_k - y_1\right) S_k \\
          \sum \left(y_k - y_1\right) T_k
      \end{bmatrix}

   This is inconsistent with the explicit note at the beginning of the
   :ref:`section <reei1-sec3>`, and with equations :eq:`gauss-pdf-S`,
   :eq:`gauss-pdf-T` and :eq:`gauss-pdf-resid`.

.. [errata-reei-4] The original expression for :math:`\sigma_1` in
   :eq:`gauss-pdf-solve` was:

   .. math::
      :label: gauss-pdf-solve-old

      \sigma_1 = -\frac{1}{B_1} \sqrt{\frac{2}{\pi}}

   The correction is due to the error in the equations :eq:`gauss-pdf-eq` and
   :eq:`gauss-pdf-int`\ [errata-reei-2]_. Notice that :math:`\mu_1` is
   unaffected by the change in factor since it is a ratio.

.. [errata-reei-5] An extension of [errata-reei-3]_ to the summary section. The
   right-most matrix had :math:`y_k - y_1` replaced with :math:`f_k - f_1`.

.. [errata-reei-6] An extension of [errata-reei-2]_ to the summary section. The
   formula for :math:`\sigma_1` has been corrected.

.. [errata-reei-7] The original paper lists equation [9]. There is no equation
   [9] in the paper, and contextually, it makes sense that the reference is in
   fact to equation [1].

.. [errata-reei-8] The equation numbers have been moved back by two. The number
   in the appendix of the original paper starts with 11. However, numbers 9 and
   10 are missing from the last numbered equation (which was
   :eq:`gauss-pdf-solve`). There are four unlabeled equations in Appendix 1,
   not 2.

.. [errata-reei-9] The table was extracted from the original figure and added
   as a separate entity. A caption was added to reflect the caption of
   :numref:`reei-gauss-pdf-data`.

.. [errata-reei-10] The figure is listed as 11 in the original paper. I am
   pretty sure that is an outright typo.

.. [errata-reei-11] The original paper has the figure number listed as 1 here,
   but should be 2.


.. include:: /link-defs.rst
