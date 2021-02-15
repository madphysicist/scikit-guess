.. _ngauss-suplement:

===========================================================
N-Dimensional Generalization of Linearized Gaussian Fitting
===========================================================

Papers by [Anthony-Granick]_ and [Wan-Wang-Wei-Li-Zhang]_ naturally focus on
the two-dimensional case of Gaussian fitting presented in image analysis. In
fact, neither paper regards axis cross-coupling terms in their representations,
instead focusing on unrotated ellipsoids.

The goal of this segment is to demonstrate a solution with the cross-coupling
term in two dimensions, and then generalize to arbitrary number of dimensions.
The mathematics of this article is very basic, but the final result is both
useful and interesting for practical applications.


.. contents::
   :local:


.. _ngauss-supplement-cross:

----------------------------
Adding a Cross-Coupling Term
----------------------------

Both papers do not go past expressing a Gaussian as

.. math::
   :label: ngauss-2d-simple

   S(x, y) = \alpha e^{-\frac{1}{2} \left( \frac{(x - x_0)^2}{\sigma_x^2} + \frac{(y - y_0)^2}{\sigma_y^2} \right)}

The complete equation, with a rotated error ellipse requires an additional
cross-coupling term :math:`\rho`:

.. math::
   :label: ngauss-2d-full

   S(x, y) = \alpha e^{-\frac{1}{2 (1 - \rho^2)} \left( \frac{(x - x_0)^2}{\sigma_x^2} - 2 \rho \frac{(x - x_0)(y - y_0)}{\sigma_x\sigma_y} + \frac{(y - y_0)^2}{\sigma_y^2} \right)}

This particular formulation is derived by inverting a 2x2 covariance matrix
with cross-coupling term :math:`\rho`:

.. math::
   :label: ngauss-2d-cov

   \Sigma_2 = \begin{bmatrix}
           \sigma_x^2             & \rho \sigma_x \sigma_y \\
           \rho \sigma_x \sigma_y & \sigma_y^2
       \end{bmatrix}

The determinant of this matrix is given by

.. math::
   :label: ngauss-2d-det

   det(\Sigma_2) = (1 - \rho^2) \sigma_x^2 \sigma_y^2

The inverse covariance matrix is therefore

.. math::
   :label: ngauss-2d-inv

   \Sigma_2^{-1} = \frac{1}{(1 - \rho^2) \sigma_x^2 \sigma_y^2} \begin{bmatrix}
           \sigma_y^2              & -\rho \sigma_x \sigma_y \\
           -\rho \sigma_x \sigma_y & \sigma_x^2
       \end{bmatrix} = \frac{1}{(1 - \rho^2)} \begin{bmatrix}
           \frac{1}{\sigma_x^2}            & -\frac{\rho}{\sigma_x \sigma_y} \\
           -\frac{\rho}{\sigma_x \sigma_y} & \frac{1}{\sigma_y^2}
       \end{bmatrix}

It should be clear how this matrix translates to :eq:`ngauss-2d-full`. If not,
see the generalization in :ref:`ngauss-supplement-ndim`, which expounds the
details more exhaustively.

If we linearize :eq:`ngauss-2d-full` by taking the logarithm of both sides and
aggregate the coefficients, we get

.. math::
   :label: ngauss-2d-log

   \text{log}(S_i) = A x_i^2 + B x_i y_i + C y_i^2 + D x_i + E y_i + F

Here :math:`S_i` is shorthand for :math:`S(x_i, y_i)`, and the parameters have been
combined as follows:

.. math::
   :label: ngauss-2d-param

   A = & -\frac{1}{2 (1 - \rho^2) \sigma_x^2} \\
   B = & \frac{\rho}{(1 - \rho^2) \sigma_x \sigma_y} \\
   C = & -\frac{1}{2 (1 - \rho^2) \sigma_y^2} \\
   D = & \frac{x_0}{(1 - \rho^2) \sigma_x^2} - \frac{\rho}{(1 - \rho^2) \sigma_x \sigma_y} \\
   E = & \frac{y_0}{(1 - \rho^2) \sigma_y^2} - \frac{\rho}{(1 - \rho^2) \sigma_x \sigma_y} \\
   F = & \text{log}(\alpha) - \frac{x_0^2}{2 (1 - \rho^2) \sigma_x^2} + \frac{\rho x_0 y_0}{(1 - \rho^2) \sigma_x \sigma_y} - \frac{y_0^2}{2 (1 - \rho^2) \sigma_y^2}

This equation is in a form suitable for a direct least-squares estimation of
the parameter values, as long as we have at least six data samples
:math:`(x_i, y_i, S_i)` available to ensure a fully determined problem. Based
on the papers that this algorithm is extending, we will want to add weights to
the solution as well. The exact nature of the weights is out of scope for this
article. Suffice it to say that we have individual weights :math:`w_i` for each
sample. Our goal is to find the projection that minimizes the error for the
following matrices:

.. math::
   :label: ngauss-2d-proj

   \begin{bmatrix}
           w_1 \cdot x_1^2 & w_1 \cdot x_1 \cdot y_1 & w_1 \cdot y_1^2 & w_1 \cdot x_1 & w_1 \cdot y_1 & w_1 \\
           w_2 \cdot x_2^2 & w_2 \cdot x_2 \cdot y_2 & w_2 \cdot y_2^2 & w_2 \cdot x_2 & w_2 \cdot y_2 & w_2 \\
           \vdots          & \vdots                  & \vdots          & \vdots        & \vdots        & \vdots \\
           w_n \cdot x_n^2 & w_n \cdot x_n \cdot y_n & w_n \cdot y_n^2 & w_n \cdot x_n & w_n \cdot y_n & w_n
       \end{bmatrix}
   \begin{bmatrix} A \\ B \\ C \\ D \\ E \\ F \end{bmatrix} =
   \begin{bmatrix}
           w_1 \cdot \text{log}(S_1) \\
           w_2 \cdot \text{log}(S_2) \\
           \vdots \\
           w_n \cdot \text{log}(S_n)
       \end{bmatrix}

Any of the common linear algebra solvers should be able to solve this least
squares problem directly.

Rather than attempting to express the six parameters :math:`x_0`, :math:`y_0`,
:math:`\sigma_x`, :math:`\sigma_y`, :math:`\rho` and :math:`A` directly in
terms of :math:`A`, :math:`B`, :math:`C`, :math:`D`, :math:`E`, :math:`F`, let
us search for a solution for :math:`x_0`, :math:`y_0`, and the elements of
:math:`\Sigma_2^{-1}` as given in :eq:`ngauss-2d-inv` instead. Not only will
this be simpler, but the generalization to mutliple dimensions will be more
apparent.

.. math::
   :label: ngauss-2d-sol1

   \Sigma_2^{-1} = -\begin{bmatrix}
           2A &  B \\
            B & 2C
       \end{bmatrix}

We can notice that the coefficients for :math:`x_0` and :math:`y_0` in the
equations for :math:`D` and :math:`E` have the form

.. math::
   :label: ngauss-2d-step2

   \begin{bmatrix} D \\ E \end{bmatrix} = & \begin{bmatrix}
           a & b \\
           b & c
       \end{bmatrix} \begin{bmatrix} x_0 \\ y_0 \end{bmatrix} \\
   a = & \frac{1}{(1 - \rho^2) \sigma_x^2} \\
   b = & -\frac{\rho}{(1 - \rho^2) \sigma_x \sigma_y} \\
   c = & \frac{1}{(1 - \rho^2) \sigma_y^2}

Inverting this matrix and working out the terms divided by the determinant
gives us

.. math::
   :label: ngauss-2d-sol2

   \begin{bmatrix} x_0 \\ y_0 \end{bmatrix} = \Sigma_2 \begin{bmatrix} D \\ E \end{bmatrix}

Finally, we can solve for the amplitude :math:`\alpha` by rewriting the
equation for :math:`F` as

.. math::
   :label: ngauss-2d-step3

   F = \text{log}(\alpha) - \frac{1}{2}\left(
           \begin{bmatrix} x_0 & y_0 \end{bmatrix}
           \Sigma_2
           \begin{bmatrix} x_0 \\ y_0 \end{bmatrix}
       \right)

The amplitude is given by

.. math::
   :label: ngauss-2d-sol3

   \alpha = e^{F + \frac{1}{2}\left(
           \begin{bmatrix} x_0 & y_0 \end{bmatrix}
           \Sigma_2
           \begin{bmatrix} x_0 \\ y_0 \end{bmatrix}
       \right)}

We can extract :math:`\sigma_x`, :math:`\sigma_y`, and :math:`\rho` from
:math:`\Sigma_2`, but as the next section shows, this is not practically
necessary.


.. _ngauss-supplement-ndim:

---------------------------------
Expanding to Arbitrary Dimensions
---------------------------------

A multivariate Gaussian is characterized by its amplitude :math:`\alpha`,
covariance matrix :math:`\Sigma`, and location :math:`\vec{\mu}`:

.. math::
   :label: ngauss-nd-full

   S(\vec{x}) = \alpha e^{-\frac{1}{2} \left( \vec{x} - \vec{\mu} \right)^T \Sigma^{-1} \left( \vec{x} - \vec{\mu} \right)}

Since :math:`\Sigma` is the positive definite matrix, we only need to specify
the upper half of it. A least squares fit to an N-dimensional Gaussian will
therefore always have :math:`\frac{N (N + 1)}{2}` parameters from the
covariance matrix, :math:`N` from the location, and one from the amplitude, for
a total of :math:`\frac{(N + 1) (N + 2)}{2}` parameters. This is consistent
with the six-parameter fit for a 2-dimensional Gaussian show in the previous
section.

We can group our coefficients into two parameter vectors and a scalar:
:math:`\vec{P}` is the vector of coefficients of the covariance terms,
:math:`\vec{Q}` is the vector of coefficients for the location terms, and
:math:`R` determines the amplitude. For the two dimensional case shown in
:eq:`ngauss-2d-param`, we have

.. math::
   :label: ngauss-nd-param

   \vec{P_2} = & \begin{bmatrix} A \\ B \\ C \end{bmatrix} \\
   \vec{Q_2} = & \begin{bmatrix} D \\ E \end{bmatrix} \\
   R_2 = & F

Since :math:`\vec{x}` is a multidimensional quantity, let us denote the
:math:`i`\ th component with a left subscript: :math:`{}_i x`. The matrix
equation that generalizes :eq:`ngauss-2d-proj` for N dimensions then becomes

.. math::
   :label: ngauss-nd-proj

   \begin{bmatrix}
           w_1 \cdot {}_1 x_1^2 & w_1 \cdot {}_1 x_1 \cdot {}_2 x_1 & \hdots &
               w_1 \cdot {}_2 x_1^2 &
               w_1 \cdot {}_2 x_1 \cdot {}_3 x_1 & \hdots &
               w_1 \cdot {}_N x_1^2 & w_1 \cdot {}_1 x_1 & \hdots &
               w_1 \\
           w_2 \cdot {}_1 x_2^2 & w_2 \cdot {}_1 x_2 \cdot {}_2 x_2 & \hdots &
               w_2 \cdot {}_2 x_2^2 &
               w_2 \cdot {}_2 x_2 \cdot {}_3 x_2 & \hdots &
               w_2 \cdot {}_N x_2^2 & w_2 \cdot {}_1 x_2 & \hdots &
               w_2 \\
           \vdots               & \vdots                            &        &
               \vdots               &
               \vdots                            &        &
               \vdots               & \vdots             &        &
               \vdots\\
           w_n \cdot {}_1 x_n^2 & w_n \cdot {}_1 x_n \cdot {}_2 x_n & \hdots &
               w_n \cdot {}_2 x_n^2 &
               w_n \cdot {}_2 x_n \cdot {}_3 x_n & \hdots &
               w_n \cdot {}_N x_n^2 & w_n \cdot {}_1 x_n & \hdots &
               w_n
       \end{bmatrix}
   \begin{bmatrix} \vec{P} \\ \vec{Q} \\ R \end{bmatrix} =
   \begin{bmatrix}
           w_1 \cdot \text{log}(S_1) \\
           w_2 \cdot \text{log}(S_2) \\
           \vdots \\
           w_n \cdot \text{log}(S_n)
       \end{bmatrix}

The solution is then a generalized form of each of :eq:`ngauss-2d-sol1`,
:eq:`ngauss-2d-sol2` and :eq:`ngauss-2d-sol3`. First we ravel :math:`\vec{P}`
to make the inverse of :math:`\Sigma`:

.. math::
   :label: ngauss-nd-sol1

   \Sigma_N^{-1} = -\begin{bmatrix}
           2 P_1  & P_2       & P_3      & \hdots & P_N \\
           P_2    & 2 P_{N+1} & P_{N+2}  & \hdots & P_{2N-1} \\
           P_3    & P_{N+2}   & 2 P_{2N} & \hdots & P_{3N-3} \\
           \vdots & \vdots    & \vdots   & \ddots & \vdots \\
           P_N    & P_{2N-1}  & P_{2N-3} & \hdots & 2 P_{\frac{N(N+1)}{2}}
       \end{bmatrix}

We can then find the location from the :math:`\vec{Q}` portion:

.. math::
   :label: ngauss-nd-sol2

   \vec{\mu} = \Sigma_N \vec{Q}

And finally the amplitude:

.. math::
   :label: ngauss-nd-sol3

   \alpha = e^{R + \frac{1}{2}\left( \vec{\mu}^T \Sigma_N \vec{\mu} \right) }
