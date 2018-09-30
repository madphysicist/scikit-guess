.. rst-class:: center

===================================
REGRESSIONS et EQUATIONS INTEGRALES
===================================

.. rst-class:: center

**Jean Jacquelin**

.. rst-class:: center

| Sample applications to various functions:
| :ref:`Gaussian <reei1-sec3>`
| :ref:`Power, Exponential, Logarithmic, Weibull <reei2>`
| :ref:`Sinusoidal <reei3>`
| Logistic
| Generalization of Sinusoidal Regression
| Damped Sinusoidal Regression
| Sum of Two Exponentials and Sum of Two Powers
| Multivariate Regression

.. image:: /_static/images/Regressions-et-equations-integrales.cover.jpg
   :name: cover-image
   :scale: 90%
   :class: center

[ First Draft : 14 January 2009 - Published : 3 January 2014 ]


.. include:: ../page_break.rst


.. rst-class:: center

-----------------
Translator's Note
-----------------

I came across this paper while searching for efficient, preferably
non-iterative, routines for fitting exponentials. The techniques presented in
this paper have served my purpose well, and this translation (and to some
extent the entire scikit) were the result. I hope that my endeavors do justice
to Jean Jacquelin's work, and prove as useful to someone else as they did to
me.

The paper translated here is a compilation of related original papers by the
author, gathered into a single multi-chapter unit. Some of the original
material is in French and some in English. I have attempted to translate the
French as faithfully as I could. I have also attempted to conform the English
portions to what I consider to be modern American usage.

A small number of technical corrections were made to the content of this paper
throughout the translation. In all cases, the errata are clearly marked with
footnotes with detailed descriptions of the fix in question. Errata are to be
submitted to Jean Jacquelin before publication of this translation.

To the extent possible, I have re-generated the plots in the paper using
`matplotlib`_ as faithfully as I could, rather than copying them out of the
original.

    -- Joseph Fox-Rabinovitz
       23rd September 2018

.. include:: /link-defs.rst

.. include:: ../page_break.rst


.. rst-class:: center

.. _reei1:

----------------------------------
Regressions and Integral Equations
----------------------------------

.. rst-class:: center

**Jean Jacquelin**


.. _reei1-abstract:

Abstract
========

The primary aim of this publication is to draw attention to a rarely used
method for solving certain types of non-linear regression problems.

The key to the method presented in this paper is the principle of linearization
through differential and integral equations, which enables the conversion of a
complicated non-linear problem into simple linear regression.

The calculus shown shown here has a fundamental difference from traditional
solutions to similar problems in that it is non-recursive, and therefore does
not require the usual iterative approach.

In order to demonstrate the theory through concrete application, detailed
numerical examples have been worked out. Regressions of power, exponential,
and logarithmic functions are presented, along with the Gaussian and Weibull
distributions commonly found in statistical applications.


.. include:: ../page_break.rst


.. rubric:: Regressions and Integral Equations
   :name: reei1-paper
   :class: center

.. rst-class:: center

**Jean Jacquelin**

.. rst-class:: center

| The first revision of the paper *Regressions and Integral Equations* was
  dated 01/14/2009.
| The current version was published on 04/27/2009.


.. _reei1-sec1:

1. Introduction
===============

The study presented here falls into the general framework of regression
problems. For example, given the coordinates of a sequence of :math:`n` points:
:math:`(x_1, y_1), (x_2, y_2), ..., (x_k, y_k), ..., (x_n, y_n)`, we wish to
find the function :math:`y = F(a, b, c, ...; x)` which lies as close as
possible to the sequence by optimizing the parameters :math:`a, b, c, ...`

The commonly known solution to linear regression merits only a brief
discussion, which is to be found in :ref:`Appendix 1 <reei1-appendix1>`. Some
problems can be solved through linear regression even though they appear
non-linear at first glance. The Gaussian distribution is an example of such a
function, and is discussed in :ref:`Appendix 2 <reei1-appendix2>`.

Excepting such simple cases, we are confronted with the daunting problem of
non-linear regression. The literature on the subject is quite extensive. Even
the briefest review would derail us from the purpose of this paper. It is also
unnecessary because our goal is to reduce some non-linear problems to linear
regression through non-iterative and non-recursive procedures (otherwise, how
would our method be innovative with respect to existing methods?).

In the next section, we will proceed to the heart of the matter: that is to
say, to render a non-linear problem to a linear form by means of a suitable
differential and/or integral equation. The preliminary discussion shows that in
the context of such problems, integral equations tend to be more numerically
stable than differential equations.

The principle of using integral equations will be explained and demonstrated in
practice using the Gaussian distribution as a concrete example. Other examples
of regression using integral equations will be described in a detailed manner
in the two following papers:

- :ref:`reei2`
- :ref:`reei3`


.. _reei1-sec2:

2. Principle of Linearization Through Differential and/or Integral Equations
============================================================================

Let us begin with a summary of numerical methods for approximating derivatives
and integrals. Given :math:`n` points :math:`(x_k, y_k)` located near the curve
of a function :math:`y(x)`, and given another function :math:`g(x)`, we can
calculate approximations for the following derivatives and integrals with
:math:`g_k = g(x_k)`:

    .. math::

       D_k = \frac{g_{k+1}y_{k+1} - g_{k-1}y_{k-1}}{x_{k+1} - x_{k-1}} \simeq
       \left(\frac{d}{dx} g(x)y(x) \right)_{\left(x = x_k \right)}

    .. math::

       DD_k = \frac{D_{k+1} - D_{k-1}}{x_{k+1} - x_{k-1}} \simeq
       \left(\frac{d^2}{dx^2} g(x)y(x)\right)_{\left(x = x_k\right)}

    And so on, for subsequent derivatives, as necessary.

    .. math::

       S_k \simeq \int_{x_1}^x g(u)y(u)du \quad
       \begin{cases}
           S_1 = 0 \hfill \text{and for $k = 2 \rightarrow n$:} & \\
           S_k = S_{k-1} +
           \frac{1}{2}(g_ky_k + g_{k-1}y_{k-1})(x_k - x_{k-1}) &
       \end{cases}

    .. math::

       SS_k \simeq \int_{x_1}^x \left(\int_{x_1}^v g(u)y(u)du\right)dv \quad
       \begin{cases}
           SS_1 = 0 \hfill \text{and for $k = 2 \rightarrow n$:} & \\
           SS_k = SS_{k-1} + \frac{1}{2}(S_k + S_{k-1})(x_k - x_{k-1}) &
       \end{cases}

    And do on, for subsequent integrals, as necessary.

It goes without saying that the points must first be ranked in order of
increasing :math:`x_k`.

It is possible to use more sophisticated approximations for numerical
differentiation and integration. Nothing prevents us from selecting a lower
limit (or lower limits) of integration other than :math:`x_1`, and using
different limits for the successive integrations. However, that would
complicate the formulas and explanations unnecessarily. For the sake of
simplicity, let us agree to use these formulas, at least for this stage of the
presentation.

Returning to our initial formulation of the problem, we wish to optimize the
parameters :math:`a, b, c, ...` of a function :math:`y(a, b, c, ...; x)` so
that its curve approaches the :math:`n` points :math:`(x_k, y_k)` as closely as
possible. Evidently, the exact expressions of the derivatives and
anti-derivatives of the function depend on the pameters :math:`a, b, c, ...`.
However the approximate values calculated using the formulas shown above, i.e.
the numerical values of :math:`D_k, DD_k, ..., S_k, SS_k, ...`, are computed
solely from the data points :math:`x_k, y_k`, **without requiring prior
knowledge of** :math:`a, b, c, ...`. This observation is the crux of the method
that is to be shown.

Let us suppose that the function :math:`y(a, b, c, ...; x)` is the solution to
a differential and/or integral equation of the form:

.. math::

   y(x) = A\Phi(x) + B\int G(x)y(x)dx + C\int\int H(x)y(x)dx^2 + ...
        + \alpha\frac{d}{dx}g(x)y(x) + \beta\frac{d^2}{dx^2}h(x)y(x) + ...

with :math:`\Phi(x), G(x), H(x), ..., g(x), h(x), ...` predetermined functions
independent of :math:`a, b, c, ...`, and
:math:`A, B, C, ..., \alpha, \beta, ...` dependent on :math:`a, b, c, ...`. The
approximate values are then respectively (with
:math:`\Phi_k = \Phi(x_k); G_k = G(x_k); H_k = H(x_k); ...`)\ [errata-reei-1]_:

.. math::

   D_k = \frac{g_{k+1}y_{k+1} - g_{k-1}y{k-1}}{x_{k+1} - x_{k-1}}

.. math::

   DD_k = \frac{\Delta_{k+1} - \Delta_{k-1}}{x_{k+1} - x_{k-1}}
   \quad \text{with} \quad
   \Delta_k = \frac{h_{k+1}y_{k+1} - h_{k-1}y_{k-1}}{x_{k+1} - x_{k-1}}

.. math::

   S_1 = 0 \quad ; \quad S_k = S_{k-1} +
   \frac{1}{2}\left(G_ky_k + G_{k-1}y_{k-1}\right)\left(x_k - x_{k-1}\right)

.. math::

   \begin{cases}
       SS_1 = 0 \quad ; \quad SS_k = SS_{k-1} +
       \frac{1}{2}\left(\Xi_k + \Xi_{k-1}\right)\left(x_k - x_{k-1}\right) \\
       \text{with: } \Xi_k = 0 \quad ; \quad \Xi_k = \Xi_{k-1} +
       \frac{1}{2}\left(H_ky_k + H_{k-1}y_{k-1}\right)\left(x_k - x_{k-1}\right)
   \end{cases}

If we replace the exact derivatives and/or anti-derivatives with their
numerical approximations, the equation will no longer hold true. We can
therefore work with the sum of squared differences:

    .. math::

       \sum_{k=1}^n \varepsilon_k^2 =
           \sum_{k=1}^n \left(-y_k + A \Phi_k + B S_k + C SS_k + ... +
           \alpha D_k + \beta DD_k + ...\right)^2

The relationship is linear with respect to
:math:`A, B, C, ..., \alpha, \beta, ...`. We have therefore returned to
classical linear regression, which allows us to calculate the optimal values of
:math:`A_0, B_0, C_0, ..., \alpha_0, \beta_0, ...`. Finally, since
:math:`A, B, C, ..., \alpha, \beta, ...` are known functions of
:math:`a, b, c, ...`, we must solve the system of equations
:math:`A(a, b, c, ...) = A_0 ; B(a, b, c, ...) = B_0 ; ... ; \alpha(a, b, c, ...) = \alpha_0 ; \beta(a, b, c, ...) = \beta_0`
to obtain the optimal values of the parameters :math:`a, b, c, ...`.

There are some additional considerations that must be taken into account when
choosing the differental and/or integral equation. Other than the requirement
for linearity in its coefficients (but not necessarily in the functions, since
we have the choice of :math:`G(x), H(x), ..., g(x), h(x), ...`), the equation
should preferably have as many coefficients
:math:`A_0, B_0, ..., \alpha_0, \beta_0, ...` as there are initial parameters
:math:`a, b, c, ...` to optimize. If there are fewer coefficients, an
additional regression (or regressions) will be necessary to calculate the
coefficients that do not figure explicitly in the equation.

Moreover, to avoid overloading the explanation, we have been considering a
simplified form of differential and/or integral equation. In fact, the equation
could have any number of different functions :math:`\Phi(x)`, several different
derivatives (corresponding to various choices of :math:`g(x)`), several
different integrals (corresponding to various choice of :math:`G(x)`), and so
on for subsequent multiple derivatives and integrals.

We see that there is a multitude of choices of differential and/or integral
equation that we can bring to bear on the problem. However, practical
considerations limit our choices. One of the main stumbling blocks is the
difficulty inherent in numerical approximation of derivatives. In fact, in
cases where the points have an irregular distribution, and are too sparse and
if, to make matters worse, the values of :math:`y_k` are not sufficiently
precise, the computed derivatives will fluctuate so much and be so dispersed as
to render the regression ineffective. On the other hand, numerical
approximations of integrals retain their stability even in these difficult
cases (this does not mean that the inevitable deviations are insignificant, but
that they remain damped, which is essential for the robustness of the overall
process). Except in special cases, the preference is therefore to seek an
integral equation rather than a differential one.

The generality of the presentation that has just been made may give the
impression that the method is complicated and difficult to implement in
practice. The fact of the matter is quite the opposite, as we will see once we
shift focus from an the abstract discussion of many possibilities to solving a
single concrete example.

One of the most spectacular examples is that of sinusoidal regression (which
we only mention in passing, without going into depth here, but which will be
treated in detail in the attached paper :ref:`reei3`). It concerns the
optimization of the parameters :math:`a, b, c` and :math:`\omega` in the
equation:

    .. math::

       y(x) = a + b sin(\omega x) + c cos(\omega x)

This function is the solution to the differential equation:

     .. math::

        y(x) = A + B \frac{d^2y}{dx^2} \quad \text{with} \quad
            A = a \quad \text{and} \quad B = -\frac{1}{\omega^2}

This is a linear equation with respect to :math:`A` and :math:`B`, which are
themselves (very simple) functions of :math:`a` and :math:`\omega`. Moreover,
the parameters :math:`b` and :math:`c` no longer figure in the differential
equation directly. This case is therefore a typical application of the method,
and among the easiest to implement, except that it contains a second
derivative, which makes it immediately suspect. Fortunately, there is no
a-priori indication that it would be better to use an integral equation whose
solution is a sinusoid instead. This method is hardly complicated and gives
largely satisfactory results (which are studied in detail in the attached
paper: :ref:`reei3`).

As a first demonstration of the calculation, let us look for a simpler example.
In the following section, we will apply the method of regression through an
integral equation to the Gaussian probability density function.


.. _reei1-sec3:

3. Example: Illustration of the Gaussian Probability Density Function
=====================================================================

We will consider a probability density function of two parameters,
:math:`\sigma` and :math:`\mu`, defined by

    .. math::
       :label: gauss-pdf-fx

       f(x) = \frac{1}{\sigma \sqrt{2 \pi}}
           exp\left(-\frac{1}{2}\left(\frac{x - \mu}{\sigma}\right)^2\right)

The general notation :math:`y(x)` of the previous sections is replaced with
:math:`f(x)` here due to the specific nature of this case.

The integration :eq:`gauss-pdf-int` leads to the integral equation
:eq:`gauss-pdf-eq`, of which :math:`f(x)` is the solution\ [errata-reei-2]_:

    .. math::
       :label: gauss-pdf-int

       \int_{x_1}^x \left(t - \mu\right)f(t)dt =
           -\sigma^2\left(f(x) - f(x_1)\right)

    .. math::
       :label: gauss-pdf-eq

       \begin{cases}
           f(x) - f(x_1) = A \int_{x_1}^x f(t)dt + B \int_{x_1}^x t f(t)dt \\
           \text{with:} \quad A = \frac{\mu}{\sigma^2} \quad
           \text{and} \quad B = -\frac{1}{\sigma^2}
       \end{cases}

This is a linear integral equation, consisting of two simple integrals, which
places it among the extensions mentioned in the previous section. We compute
the respective approximations, the first denoted :math:`S` with
:math:`G(x) = 1`, and the second denoted :math:`T` with :math:`G(x) = x`:

    .. math::
       :label: gauss-pdf-S

       \begin{cases}
           S_1 = 0 \\
           S_k = S_{k-1} +
               \frac{1}{2}\left(f_k + f_{k-1}\right)
               \left(x_k - x_{k-1}\right) \quad k = 2 \rightarrow n
       \end{cases}

    .. math::
       :label: gauss-pdf-T

       \begin{cases}
           T_1 = 0 \\
           T_k = T_{k-1} +
               \frac{1}{2}\left(x_k f_k + x_{k-1} f_{k-1}\right)
               \left(x_k - x_{k-1}\right) \quad k = 2 \rightarrow n
       \end{cases}

When we replace :math:`f(x_k)` with :math:`f_k`, :math:`f(x_1)` with
:math:`f_1` and the integrals with :math:`S_k` and :math:`T_k`, respectively,
equation :eq:`gauss-pdf-eq` no longer holds true. We seek to minimize the sum of
the squares of the residuals:

    .. math::
       :label: gauss-pdf-resid

       \sum_{k=1}^n \varepsilon_k^2 =
           \sum_{k=1}^n \left(-\left(f_k - f_1\right) + A S_k + B T_k\right)^2

Notice that if we had chosen a different lower limit for the integration than
:math:`x_1`, it would have changed not only the value of :math:`f_1`, but also
the numerical values of :math:`S_k` and :math:`T_k`, in a way that would cancel
out the differences without changing the final result.

The relationship :eq:`gauss-pdf-resid` is none other than the than the equation
of a linear regression, which we know how to optimize for the parameters
:math:`A_1` and :math:`B_1`\ [errata-reei-3]_:

    .. math::
       :label: gauss-pdf-lsq

       \begin{bmatrix}A_1 \\ B_1\end{bmatrix} =
       \begin{bmatrix}
           \sum \left(S_k\right)^2 & \sum S_k T_k \\
           \sum S_k T_k            & \sum \left(T_k\right)^2
       \end{bmatrix}
       \begin{bmatrix}
           \sum \left(f_k - f_1\right) S_k \\
           \sum \left(f_k - f_1\right) T_k
       \end{bmatrix}

With the convention that :math:`\sum \equiv \sum_{k=1}^n`. We then deduce
:math:`\sigma_1` and :math:`\mu_1` according to :eq:`gauss-pdf-eq`\ [errata-reei-4]_:

     .. math::
        :label: gauss-pdf-solve

        \sigma_1 = \sqrt{-\frac{1}{B_1}} \quad ; \quad \mu_1 = -\frac{A_1}{B_1}

Here is a summary of the numerical computation:

.. _reei1-sec3-alg:

.. rst-class:: listing

.. 

   **Data:** :math:`(x_1, f_1), (x_2, f_2), ..., (x_k, f_k), ..., (x_n, f_n)`

   - Compute :math:`S_k`:

     .. math::

        \begin{cases}
            S_1 = 0 \\
            S_k = S_{k-1} +
                \frac{1}{2}\left(f_k + f_{k-1}\right)
                \left(x_k - x_{k-1}\right) \quad k = 2 \rightarrow n
        \end{cases}

   - Compute :math:`T_k`:

     .. math::

        \begin{cases}
            T_1 = 0 \\
            T_k = T_{k-1} +
                \frac{1}{2}\left(x_k f_k + x_{k-1} f_{k-1}\right)
                \left(x_k - x_{k-1}\right) \quad k = 2 \rightarrow n
        \end{cases}

   - Compute :math:`\sum \left(S_k\right)^2 , \sum S_k T_k, \sum \left(T_k\right)^2 \\ \sum \left(f_k - f_1\right) S_k, \sum \left(f_k - f_1\right) T_k`

   - Compute :math:`A_1` and :math:`B_1`\ [errata-reei-5]_:

     .. math::

        \begin{bmatrix}A_1 \\ B_1\end{bmatrix} =
        \begin{bmatrix}
            \sum \left(S_k\right)^2 & \sum S_k T_k \\
            \sum S_k T_k            & \sum \left(T_k\right)^2
        \end{bmatrix}
        \begin{bmatrix}
            \sum \left(f_k - f_1\right) S_k \\
            \sum \left(f_k - f_1\right) T_k
        \end{bmatrix}

   - Compute :math:`\sigma_1` and :math:`\mu_1`\ [errata-reei-6]_:
     :math:`\sigma_1 = \sqrt{-\frac{1}{B_1}} \quad ; \quad \mu_1 = -\frac{A_1}{B_1}`

   **Result:** :math:`\sigma_1` and :math:`\mu_1` are the approximations of
   :math:`\sigma` and :math:`\mu`.

To illustrate the calculation (:numref:`reei-gauss-pdf-plot`), numerical data
(:numref:`reei-gauss-pdf-data`) was genererated in the following manner:
:math:`x_k` values were chosen at random from the domain. From the "exact"
values :math:`\sigma_e` and :math:`\mu_e` (defining the "exact" function
:math:`f(x)`, whose representative curve is plotted as a dashed line in
:numref:`reei-gauss-pdf-plot`), we computed the exact values of :math:`f(x_k)`
given by equation :eq:`gauss-pdf-fx`\ [errata-reei-7]_. Then we added
deviations randomly drawn from a range - to +10% of :math:`f(x_k)`, which gave
us the numerical values of :math:`f_k` in :numref:`reei-gauss-pdf-data`, after
rounding.

The outrageous error modeling is motivated by the need for legibility in the
figure, so that the so called "experimental" points, represented by crosses,
lie far enough from the dashed curve to be clearly distinguishable. In the same
vein, only a handful of points was chosen so that the deviations between the
"exact" dashed curves and the calculated solid curves are are highlighted for
both the intermediate and the final calculation. The fact that the points are
not uniformly distributed across the domain is also a significant complication.

.. figure:: /generated/reei/gauss-pdf-plot.png
   :name: reei-gauss-pdf-plot

   A sample regression of the Gaussian probability density function.

.. table:: Numerical values corresponding to the example in \
   :numref:`reei-gauss-pdf-plot`.
   :name: reei-gauss-pdf-data
   :class: data-table

   .. include:: /generated/reei/gauss-pdf-data.rst

In :numref:`reei-gauss-pdf-plot`, the shape of the curves of the "exact"
integrals and the points :math:`\left(x_k, S_k\right)` and
:math:`\left(x_k, T_k\right)` make the primary reason for the deviations in
this method of calculation clearly apparent: numerical integration, while
preferable to derivation, is not by any means perfect, and causes the
deviations in :math:`\left(\sigma_1, \mu_1\right)`.

To form an objective opinion about the qualities and defects of the method
exposed here, it would be necessary to perform a systematic experimental
study of a very large number of cases and examples. In the current state of
progress of this study, this remains yet to be done.

It is certain that the deviations, caused by the defects inherent in numerical
integration, would be considerably reduced if the points were sufficiently
numerous and their abscissae were partitioned at sufficiently regular
intervals.


.. _reei1-sec4:

4. Comments
===========

It would be unreasonable to imagine that the method presented here could
replace currently used methods, implemented in commercial software, and with
the benefits of long term study, experimentation and improvements. We must ask
then, what is the motivation behind this work.

Of course, recursive methods generally require a first approximation within the
same order of magnitude as the target value. This is not generally a handicap
since users are not completely in the dark. One might be tempted to consider
this method of integral equations to satisfy the need for initial
approximation. However, the need is quite marginal, and should not be seen as a
serious motivation.

A simple method, easy to program, like the one shown here, might certainly
seduce some potential users in particular situations where we seek full mastery
over the calculations that are performed. Users of commercial software are
usually satisfied with the results they provide, but may occasionally deplore
not knowing precisely what the sophisticated software is doing. Nevertheless,
it would be poor motivation indeed for this study to attempt to provide a less
powerful tool than what already exists, just to aleviate some of the
frustration caused by using tools whose exact mechanisms are unknown.

In fact, we must see in this paper not a specific utilitarian motivation in
the case of the Gaussian distribution, but the intention to understand and draw
attention to a more general idea: the numerous possibilities offered by
integral equations to transform a problem of non-linear regression into a
linear regression, computed through a non-iterative process.

It is out of the question to compete against what already exists, and what is
more imporant, works well. On the other hand, it would be a pity to forget a
method that might potentially help resolve future problems: the method that has
been the subject of this paper, whose essence is presented in
:ref:`section 2 <reei1-sec2>`.


.. rst-class:: center

.. _reei1-appendix1:

Appendix 1: Review of Linear Regression
=======================================

When the function that we seek to optimize, :math:`y = F(a, b, c, ...; x)`, can
be written in the form :math:`y = a f(x) + b g(x) + c h(x) + ...`, according to
the number of parameters :math:`a, b, c, ...` and with the functions
:math:`f(x), g(x), h(x), ...` being known, the process is linear with respect
to the optimization parameters.

Even more generally, if the function :math:`y = F(a, b, c, ...; x)` can be
transformed into :math:`F(x, y) = A f(x, y) + B g(x, y) + C h(x, y)` with known
functions
:math:`F(x, y), f(x, y), g(x, y), h(x, y), ..., A(a, b, c, ...), B(a, b, c, ...), C(a, b, c, ...), ...`
the process is again linear with respect to the coefficients :math:`A`,
:math:`B` and :math:`C`, even if it is not linear with respect to
:math:`a, b, c, ...`. But it always reverts to a linear regression. The method
of "least squares" effectively consists of finding the minimum of:

    .. math::

       \begin{cases}
           \varepsilon^2_{\left(A, B, C, ...\right)} =
               \sum_{k=1}^n \left(F_k -
                   \left(A f_k + B g_k + C h_k + ...\right)\right)^2 \\
           F_k \equiv F(x_k, y_k) ; f_k \equiv f(x_k, y_k) ;
               g_k \equiv g(x_k, y_k) ; h_k \equiv h(x_k, y_k)
       \end{cases}

The partial derivatives with respect to :math:`A, B, C, ...` determine a system
of equations whose solutions, :math:`A_0, B_0, C_0, ...` are optimal:

    .. math::

       \begin{cases}
           \left(\frac{\partial \left(\varepsilon^2\right)}{\partial A}\right)
                   _{A_0, B_0, C_0, ...} = -\sum_{k=1}^n \left(F_k - \left(
               A_0 f_k + B_0 g_k + C_0 h_k, ...\right)\right) f_k = 0 \\
           \left(\frac{\partial \left(\varepsilon^2\right)}{\partial B}\right)
                   _{A_0, B_0, C_0, ...} = -\sum_{k=1}^n \left(F_k - \left(
               A_0 f_k + B_0 g_k + C_0 h_k, ...\right)\right) g_k = 0 \\
           \left(\frac{\partial \left(\varepsilon^2\right)}{\partial C}\right)
                   _{A_0, B_0, C_0, ...} = -\sum_{k=1}^n \left(F_k - \left(
               A_0 f_k + B_0 g_k + C_0 h_k, ...\right)\right) h_k = 0 \\
           ...
       \end{cases}

The solution to this system, conventionally written with
:math:`\sum \equiv \sum_{k=1}^n` leads to:

    .. math::

       \begin{bmatrix}
           A_0 \\ B_0 \\ C_0 \\ ...
       \end{bmatrix} =
       \begin{bmatrix}
           \sum f_k^2   & \sum f_k g_k & \sum f_k h_k & ... \\
           \sum f_k g_k & \sum g_k^2   & \sum g_k h_k & ... \\
           \sum f_k h_k & \sum g_k h_k & \sum h_k^2   & ... \\
           ...          & ...          & ...          & ...
       \end{bmatrix}^{-1}
       \begin{bmatrix}
           \sum F_k f_k \\ \sum F_k g_k \\ \sum F_k h_k \\ ...
       \end{bmatrix}

Then we obtain the optimized values of :math:`a, b, c, ...` corresponding to
the following system, where the unknowns are :math:`a_0, b_0, c_0, ...`:

    .. math::

       \begin{cases}
           A(a_0, b_0, c_0, ...) = A_0 \\
           B(a_0, b_0, c_0, ...) = B_0 \\
           C(a_0, b_0, c_0, ...) = C_0 \\
           ...
       \end{cases}

which is a system that is non-linear in the same measure that the functions
:math:`A(a, b, c, ...), B(a, b, c, ...), C(a, b, c, ...), ...` are non-linear.
But this does not prevent the regression that was performed from being linear,
so even this case has its rightful place in this section.

Of course, this can be further extended by considering more variables, for
example :math:`x, y, z, t, ...`, instead of just :math:`x, y`, thus working in
3D, or 4D, ..., instead of 2D. Everything mentioned here figures in literature
in a more detailed and more structured manner, with presentations adapted to
the exposition of general theory. The purpose was only to present a brief
review, with a specific notation to be used consistently throughout the
remainder of the work.


.. rst-class:: center

.. _reei1-appendix2:

Appendix 2: Linear Regression of the Gaussian Cumulative Distribution Function
==============================================================================

We consider the cumulative Gaussian distribution function of two parameters,
:math:`\sigma` and :math:`\mu`, defined by\ [errata-reei-8]_:

    .. math::
       :label: gauss-cdf-fx

       F(x) = \frac{1}{\sqrt{2 \pi} \sigma} \int_{-\inf}^x exp\left(
           -\frac{1}{2}\left(\frac{t - \mu}{\sigma}\right)^2\right)dt

An example is shown in :numref:`reei-gauss-cdf-plot` (the dashed curve).

.. figure:: /generated/reei/gauss-cdf-plot.png
   :name: reei-gauss-cdf-plot

   Sample regression illustrating the Gaussian cumulative distribution function.

.. table:: Numerical values corresponding to the example in \
   :numref:`reei-gauss-cdf-plot`\ [errata-reei-9]_.
   :name: reei-gauss-cdf-data
   :class: data-table

   .. include:: /generated/reei/gauss-cdf-data.rst

The data are the points that we call "experimental":
:math:`(x_1, F_1), (x_2, F_2), ..., (x_k, F_k), ..., (x_n, F_n)`, which, in
:numref:`reei-gauss-cdf-plot`\ [errata-reei-10]_ have a particular dispersion
relative to their respective nominal positions :math:`\left(x_k, F(x_k)\right)`
on the dashed curve representing :math:`F(x)`.


.. _reei1-appendix2-listings:

Listings for the functions Erf and argErf
-----------------------------------------


.. rst-class:: center

.. _reei2:

------------------------------------------------------------------------------
Non-Linear Regression of Power, Exponential, Logarithmic and Weibull Functions
------------------------------------------------------------------------------

.. rst-class:: center

**Jean Jacquelin**


.. _reei2-abstract:

Abstract
========

We demonstrate the application of a well-chosen integral equation to produce a
non-iterative optimization of the parameters of power, exponential, logarithmic
and Weibull functions.


.. rubric:: Non-Linear Regression of Power, Exponential, Logarithmic and
            Weibull Functions
   :name: reei2-paper
   :class: center


.. rst-class:: center

**Jean Jacquelin**


.. _reei2-sec1:

1. Introduction
===============


.. _reei2-sec2:

2. Regression of Functions of the Form :math:`y(x) = a + b exp(c x)`
====================================================================


.. _reei2-sec3:

3. Regression of the Three-Parameter Weibull Distribution
=========================================================


.. _reei2-sec4:

4. Conclusion
=============


.. rst-class:: center

.. _reei3:

-----------------------
Regression of Sinusoids
-----------------------

.. rst-class:: center

**Jean Jacquelin**

.. rst-class:: center

| The first revision of the paper *Regressions of Sinusoids* was dated
  01/09/2009.
| The current version was published on 02/15/2009.


.. _reei3-sec1:

1. Introduction
===============


.. _reei3-sec2:

2. Instances Where :math:`\omega` is Known Ahead of Time
========================================================


.. _reei3-sec3:

3. Linearization With an Integral Equation
==========================================


.. _reei3-sec4:

4. Succinct Performance Analysis
================================


.. _reei3-sec4-1:

4.1 "Equidistant" Distribution of Abscissae and Non-dispersion of Ordinals
--------------------------------------------------------------------------


.. _reei3-sec4-2:

4.2 Aleatory Distribution of Point Abscissae Without Ordinal Dispersion
-----------------------------------------------------------------------


.. _reei3-sec4-3:

4.3 Aleatory Distribution of Point Abscissae With Dispersed Ordinals
--------------------------------------------------------------------


.. _reei3-sec5:

5. Instances Where :math:`a` and :math:`\rho` Parameters Are Approximately Known
================================================================================


.. _reei3-sec6:

6. Results of a Complete Optimization
=====================================


.. _reei3-sec7:

7. Commentary
=============


.. rst-class:: center

.. _reei3-appendix1:

Appendix 1: Summary of Sinusoidal Regression Algorithm
======================================================


.. rst-class:: center

.. _reei3-appendix2:

Appendix 2: Detailed Procedure for Sinusoidal Regression
========================================================


.. rst-class:: center

.. _reei4:

-----------------------------------------------------------
Application to the Logistic Distribution (Three Parameters)
-----------------------------------------------------------



.. rst-class:: center

.. _reei5:

----------------------------------------------------------
Application to the Logistic Distribution (Four Parameters)
----------------------------------------------------------


.. rst-class:: center

.. _reei6:

--------------------------------------
Mixed Linear and Sinusoidal Regression
--------------------------------------


.. rst-class:: center

.. _reei7:

---------------------------------
Generalized Sinusoidal Regression
---------------------------------


.. rst-class:: center

.. _reei8:

----------------------------
Damped Sinusoidal Regression
----------------------------


.. rst-class:: center

.. _reei9:

-------------------------------------------------------
Double Exponential Regression & Double Power Regression
-------------------------------------------------------


.. rst-class:: center

.. _reei10:

-----------------------
Multivariate Regression
-----------------------
