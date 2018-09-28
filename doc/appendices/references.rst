==========
References
==========

This section of the documentation is devoted to providing references for the
algorithms implemented in scikit-guess. Each paper comes with a link, a PDF
where permitted, and any additional materials.

.. contents:: List of References
   :local:
..   :depth: 1


.. _ref-reei:

-----------------------------------
Régressions et équations intégrales
-----------------------------------

The ideas proposed in this paper by Jean Jacquelin were the seed for this
scikit. A PDF is available with this documentation:
:download:`Régressions et équations intégrales </_static/documents/Regressions-et-equations-integrales.pdf>`.

Available online at
https://www.scribd.com/doc/14674814/Regressions-et-equations-integrales.


Concept
=======

The concept behind this paper is that integrals and derivatives can be
estimated through differentials and cumulative sums. The goal is to set up an
integral or differential equation whose solution is the model function. If the
right-hand side of such an equation can be expressed as a linear combination of
integrals and derivatives of itself multiplied by some other pre-determined
functions, we can estimate the numerical values of the equation's terms. While
the coefficients that make the equation work depend on the fitting parameters,
the functions themselves do not. It is therefore possible to set up a simple
linear regression for the coefficients based on the numerical approximations of
the integrals and derivatives. The numerical approximation of integrals by
cumulative sums tend to be more robust than the approximations of derivatives
by differentials, so integral equations are generally preferred.


Translation
===========

The original paper is mostly in French, so an English translation is provided
as part of the documentation of scikit-guess. The translation can be read here:

.. toctree::
   :maxdepth: 4
   :caption: Paper Contents

   reei/translation
   reei/supplement
