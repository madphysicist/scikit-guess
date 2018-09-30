=======================
Supplementary Materials
=======================

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

   x = \left( A^T * A \right)^{-1} \left( A^T * b \right)

The formulae in the paper show how to compute the elements of :math:`A^T * A`
and :math:`A^T * b`, which is usually done more efficiently by exiting packages
when given the raw :math:`A` and :math:`b`.

The following sections show how to construct such simplified solutions to the
equations in the paper. Solutions are described conceptually, and presented
concretely with Python code. The solutions here are for conceptual reference
only. They are not a complete or exact reflection of how things are done in the
scikit itself.

In the code snippets below, ``x`` and ``y`` are assumed to be one-dimensional
numpy arrays. Both arrays have an equal number of elements, ``n``. Numpy is
implicitly imported under the conventional name ``np``.

Gaussian
========

Algorithm originally presented in :ref:`reei1-sec3` and summarized
:ref:`here <reei1-sec3-alg>`.

:math:`A` is a matrix with the cumulative sums :math:`S` and :math:`T` as
the columns. In numpy terms:

.. code-block:: python

   S = np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(x))
   S = np.insert(S, 0, 0)

   xy = x * y
   T = np.cumsum(0.5 * (xy[1:] + xy[:-1]) * np.diff(x))
   T = np.insert(T, 0, 0)

   A = np.stack((S, T), axis=1)

:math:`b` is the vector of measured :math:`y` values decreased by its first
element. In numpy terms:

.. code-block:: python

   b = y - y[0]


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


.. include:: /link-defs.rst
