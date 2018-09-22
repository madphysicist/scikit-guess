.. rst-class:: center

===================================
REGRESSIONS et EQUATIONS INTEGRALES
===================================

.. rst-class:: center

**Jean Jacquelin**

.. rst-class:: center

| Sample applications to various functions:
| :ref:`Gaussian <x1-sec3>`
| Power, Exponential, Logarithmic, Weibull
| Sinusoidal
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

The paper translated here is a compilation of eight original papers by the
author, gathered into a single multi-chapter unit. Some of the original
matrial is in French and some in English. I have attempted to translate the
French as faithfully as I could. I have also attempted to clean up the grammar
and usage of the English portions. My goal is to represent the meaning of the
original as accurately possible, both technically and linguistically. I have
made shameless use of Google Translate, some dictionaries that I have held on
to since before high school, and some lingering memories of my French classes.
In all cases, the final wording is entirely my own.


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


.. rst-class:: center

.. _x1-paper:

Regressions and Integral Equations
==================================

.. rst-class:: center

**Jean Jacquelin**

.. rst-class:: center

| The first revision of the paper *Regressions and Integral Equations* was
  dated 01/14/2009.
| The current version was published on 04/27/09.


.. _x1-sec1:

1. Introduction
---------------

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


.. _x1-sec2:

2. Principle of Linearization Through Differential and/or Integral Equations
----------------------------------------------------------------------------


.. _x1-sec3:

3. Example: Case of the Gaussian Probability Density Function
-------------------------------------------------------------


.. _x1-sec4:

4. Commentary
-------------


.. rst-class:: center

.. _x1-appendix1:

Appendix 1: Linear Regression (Refresher)
-----------------------------------------


.. rst-class:: center

.. _x1-appendix2:

Appendix 2: Case of the Gaussian Probability Density Function
-------------------------------------------------------------

