.. rst-class:: center

===================================
REGRESSIONS et EQUATIONS INTEGRALES
===================================

.. rst-class:: center

**Jean Jacquelin**

.. rst-class:: center

| Sample applications to various functions:
| :ref:`Gaussian <x1-sec3>`
| :ref:`Power, Exponential, Logarithmic, Weibull <x2>`
| :ref:`Sinusoidal <x3>`
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


.. include:: page_break.rst


.. rst-class:: center

-----------------
Translator's Note
-----------------

I came across this paper while searching for efficient optimization routines
for fitting exponentials. The techniques presented in this paper have served my
purpose well, and this translation (as well as the whole scikit) were the
result. I hope that these endeavors do justice to Jean Jacquelin's work, and
prove as useful to someone else as they did to me.

The paper translated here is a compilation of related original papers by the
author, gathered into a single multi-chapter unit. Some of the original matrial
is in French and some in English. I have attempted to translate the French as
faithfully as I could. I have also attempted to clean up the grammar and usage
of the English portions. My goal is to represent the meaning of the original as
accurately possible, both technically and linguistically.


[ Translation : 22 September 2018 ]


.. include:: page_break.rst


.. rst-class:: center

.. _x1:

----------------------------------
Regressions and Integral Equations
----------------------------------

.. rst-class:: center

**Jean Jacquelin**


.. _x1-abstract:

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


.. include:: page_break.rst


.. rubric:: Regressions and Integral Equations
   :name: x1-paper
   :class: center

.. rst-class:: center

**Jean Jacquelin**

.. rst-class:: center

| The first revision of the paper *Regressions and Integral Equations* was
  dated 01/14/2009.
| The current version was published on 04/27/2009.


.. _x1-sec1:

1. Introduction
===============

The survey presented here falls into the general category of regression
problems. For example, given the coordinates of a sequence of :math:`n` points:
:math:`(x_1, y_1), (x_2, y_2), ..., (x_k, y_k), ..., (x_n, y_n)`, we wish to
find the function :math:`y = F(a, b, c, ...; x)` which lies as close as
possible to the sequence by optimizing the parameters :math:`a, b, c, ...`

The commonly known solution to linear regression merits only a brief
discussion, which is to be found in :ref:`x1-appendix1`. There are a few cases
whose linearity might be overlooked at first glance, but for which it is
possible to recast the problem as a linear regression. The case of the Gaussian
distribution function is an example, treated in :ref:`x1-appendix2`.

Barring such trivial cases, we are confronted with the daunting problem of
non-linear regression. The literature on the subject is quite extensive. Even
the briefest review would derail us from the purpose of this paper. It is also
unnecessary because our goal is to reduce some non-linear problems to a linear
regression without using iterative or recursive procedures (otherwise, how
would our method be innovative with respect to existing methods?).

Starting with the next paragraph, we will proceed to the heart of the matter:
that is to say, to render a non-linear problem to a linear form using a
suitable differential and/or integral equation. The preliminary discussion
shows that, except in special cases, an integral equation is better suited to
the task of numerical approximation than a differential equation, in the
context of such problems.

The idea of using integral equations will be explained and demonstrated in
practice using a Gaussian distribution as a concrete example. Other examples of
regression using integral equations will be described in a detailed manner in
the two associated papers:

- :ref:`x2`
- :ref:`x3`


.. _x1-sec2:

2. Principle of Linearization Through Differential and/or Integral Equations
============================================================================

Let us begin with a summary of the numerical methods used to approximate
derivatives and/or primitives. Given :math:`n` points :math:`(x_k, y_k)` that
lie near the curve of a function :math:`y(x)`, and given another function
:math:`g(x)`, we can calculate the approximations for the following derivatives
and/or integrals, with :math:`g_k = g(x_k)`:

.. math::

   D_k = \frac{g_{k+1}y_{k+1} - g_{k-1}y_{k-1}}{x_{k+1} - x_{k-1}} \simeq
   \left(\frac{d}{dx} g(x)y(x) \right)_{\left(x = x_k \right)}

.. math::

   DD_k = \frac{D_{k+1} - D_{k-1}}{x_{k+1} - x_{k-1}} \simeq
   \left(\frac{d^2}{dx^2} g(x)y(x)\right)_{\left(x = x_k\right)}

And so on, for further derivatives, as necessary.

.. math::

   S_k \simeq \int_{x_1}^x g(u)y(u)du
   \begin{cases}
       S_1 = 0 \hfill \text{and for $k = 2 \rightarrow n$:} & \\
       S_k = S_{k-1} + \frac{1}{2}(g_ky_k + g_{k-1}y_{k-1})(x_k - x_{k-1}) &
   \end{cases}

.. math::

   SS_k \simeq \int_{x_1}^x \left(\int_{x_1}^v g(u)y(u)du\right)dv
   \begin{cases}
       SS_1 = 0 \hfill \text{and for $k = 2 \rightarrow n$:} & \\
       SS_k = SS_{k-1} + \frac{1}{2}(S_k + S_{k-1})(x_k - x_{k-1}) &
   \end{cases}

And do on, for further integrals, as necessary.

It goes without saying that the points must first be ranked in order of
increasing :math:`x_k`.


.. _x1-sec3:

3. Example: Case of the Gaussian Probability Density Function
=============================================================


.. _x1-sec4:

4. Commentary
=============


.. rst-class:: center

.. _x1-appendix1:

Appendix 1: Review of Linear Regression
=======================================


.. rst-class:: center

.. _x1-appendix2:

Appendix 2: Linear Regression of the Gaussian Probability Density Function
==========================================================================


.. rst-class:: center

.. _x2:

------------------------------------------------------------------------------
Non-Linear Regression of Power, Exponential, Logarithmic and Weibull Functions
------------------------------------------------------------------------------

.. rst-class:: center

**Jean Jacquelin**


.. _x2-abstract:

Abstract
========

We demonstrate the application of a well-chosen integral equation to produce a
non-iterative optimization of the parameters of power, exponential, logarithmic
and Weibull functions.


.. rubric:: Non-Linear Regression of Power, Exponential, Logarithmic and
            Weibull Functions
   :name: x2-paper
   :class: center


.. rst-class:: center

**Jean Jacquelin**


.. _x2-sec1:

1. Introduction
===============


.. _x2-sec2:

2. Regression of Functions of the Form :math:`y(x) = a + b exp(c x)`
====================================================================


.. _x2-sec3:

3. Regression of the Three-Parameter Weibull Distribution
=========================================================


.. _x2-sec4:

4. Conclusion
=============


.. rst-class:: center

.. _x3:

-----------------------
Regression of Sinusoids
-----------------------

.. rst-class:: center

**Jean Jacquelin**

.. rst-class:: center

| The first revision of the paper *Regressions of Sinusoids* was dated
  01/09/2009.
| The current version was published on 02/15/2009.


.. _x3-sec1:

1. Introduction
===============


.. _x3-sec2:

2. Cases Where :math:`\omega` is Known Ahead of Time
====================================================


.. _x3-sec3:

3. Linearization With an Integral Equation
==========================================


.. _x3-sec4:

4. Succinct Performance Analysis
================================


.. _x3-sec4-1:

4.1 "Equidistant" Distribution of Abscissae and Non-dispersion of Ordinals
--------------------------------------------------------------------------


.. _x3-sec4-2:

4.2 Aleatory Distribution of Point Abscissae Without Ordinal Dispersion
-----------------------------------------------------------------------


.. _x3-sec4-3:

4.3 Aleatory Distribution of Point Abscissae With Dispersed Ordinals
--------------------------------------------------------------------


.. _x3-sec5:

5. Cases Where :math:`a` and :math:`\rho` Parameters Are Approximately Known
============================================================================


.. _x3-sec6:

6. Results of a Complete Optimization
=====================================


.. _x3-sec7:

7. Commentary
=============


.. rst-class:: center

.. _x3-appendix1:

Appendix 1: Summary of Sinusoidal Regression Algorithm
======================================================


.. rst-class:: center

.. _x3-appendix2:

Appendix 2: Detailed Procedure for Sinusoidal Regression
========================================================


.. rst-class:: center

.. _x4:

-----------------------------------------------------------
Application to the Logistic Distribution (Three Parameters)
-----------------------------------------------------------



.. rst-class:: center

.. _x5:

----------------------------------------------------------
Application to the Logistic Distribution (Four Parameters)
----------------------------------------------------------


.. rst-class:: center

.. _x6:

--------------------------------------
Mixed Linear and Sinusoidal Regression
--------------------------------------


.. rst-class:: center

.. _x7:

---------------------------------
Generalized Sinusoidal Regression
---------------------------------


.. rst-class:: center

.. _x8:

----------------------------
Damped Sinusoidal Regression
----------------------------


.. rst-class:: center

.. _x9:

-------------------------------------------------------
Double Exponential Regression & Double Power Regression
-------------------------------------------------------


.. rst-class:: center

.. _x10:

-----------------------
Multivariate Regression
-----------------------
