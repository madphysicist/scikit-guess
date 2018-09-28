=======================
Supplementary Materials
=======================

------
Errata
------

The following list of errata was complied during the translation. All items
link back to the corresponding location in the translated paper as footnotes.

.. [err-x1-sec2-1] The following was omitted from the list, because I believe
   that it is not correct/necessary, pending confirmation from the author:
   :math:`; \alpha_k = \alpha(x_k); \beta_k = \beta(x_k); ...`.

.. [err-x1-sec3-1] The original equations :eq:`gauss-int` and :eq:`gauss-eq`
   read:

   .. math::
      :label: gauss-int-old

      \int_{x_1}^x \left(t - \mu\right)f(t)dt =
          -\sqrt{\frac{\pi}{2}}\sigma\left(f(x) - f(x_1)\right)

   .. math::
      :label: gauss-eq-old

      \begin{cases}
          f(x) - f(x_1) = A \int_{x_1}^x f(t)dt + B \int_{x_1}^x t f(t)dt \\
          \text{with:} \quad A = \frac{\mu}{\sigma}\sqrt{\frac{2}{\pi}} \quad
          \text{and} \quad B = -\frac{1}{\sigma}\sqrt{\frac{2}{\pi}}
      \end{cases}

   Taking the integral :math:`\int_{x_1}^x \left(t - \mu\right) f(t)dt`, with
   :math:`f(x)` given in :eq:`gauss-fx`, shows that the factor of
   :math:`\sqrt{\frac{\pi}{2}}` should actually be an additional factor of
   :math:`\sigma`. This factor then propagates to :math:`A` and :math:`B` as
   shown in the translated paper.

.. [err-x1-sec3-2] The original equation :eq:`gauss-lsq` had :math:`y - y_1`
   in the rightmost matrix:

   .. math::
      :label: gauss-lsq-old

      \begin{bmatrix}
          \sum \left(y_k - y_1\right) S_k \\
          \sum \left(y_k - y_1\right) T_k
      \end{bmatrix}

   This is inconsistent with the explicit note at the beginning of the
   :ref:`section <x1-sec3>`, and with equations :eq:`gauss-S`, :eq:`gauss-T`
   and :eq:`gauss-resid`.

.. [err-x1-sec3-3] The original expression for :math:`\sigma_1` in
   :eq:`gauss-solve` was:

   .. math::
      :label: gauss-solve-old

      \sigma_1 = -\frac{1}{B_1} \sqrt{\frac{2}{\pi}}

   The correction is due to the error in the equations :eq:`gauss-eq` and
   :eq:`gauss-int`. Notice that :math:`\mu_1` is unaffected by the change in
   factor since it is a ratio.


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

Algorithm originally presented in :ref:`x1-sec3` and summarized
:ref:`here <x1-sec3-alg>`.

:math:`A` is a concatenation of the cumulative sums :math:`S` and :math:`T` as
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


.. include:: /link-defs.rst