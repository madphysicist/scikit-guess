==========
References
==========

This section of the documentation is devoted to providing references for the
algorithms implemented in scikit-guess. Each paper comes with a link, a PDF
where permitted, and any additional materials.

.. contents:: List of References
   :local:
..   :depth: 1


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

The concept behind this paper is that any function can be rewritten as a linear
combination of integrals and derivatives of itself multiplied by some other
function. By finding the appropriate integral or differential equation, we can
estimate the parameters of the function that would optmimze its fit to our
data. We can estimate the values of the integrals and derivatives from our data
using either simple summation and differencing (as in the paper), or more
complex techniques, thereby setting up a simple system of equations for the
optimal parameters amounting to a simple least squares regression. The
approximations of integrals by cumulative sums tend to be more robust than the
approximations of derivatives by differentials, so integral equations are
preferred.


Translation
===========

The original paper is mostly in French, so an English translation is provided
as part of the documentation of scikit-guess. The translation can be read here:

.. toctree::
   :maxdepth: 4
   :caption: Paper Contents

   translation