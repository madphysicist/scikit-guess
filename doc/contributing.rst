============
Contributing
============

If you have an idea or would like to correct something, please submit an
`issue <https://github.com/madphysicist/scikit-guess/issues>`_, or a
`pull request <https://github.com/madphysicist/scikit-guess/pulls>`_.

If you have questions or are not sure if something is a bug, feel free to
submit an issue, or ask a question tagged
`scikit-guess <https://stackoverflow.com/questions/tagged/scikit-guess>`_ on
`Stack Overflow <https://stackoverflow.com/>`_.

This scikit is still in early stages, so please provide all the criticism and
advice you can. Any support at all is welcome. In particular, the following
areas would be appreciated:

- Proper testing setup (e.g. TravisCI, Appveyor)
- Pointing out anything that is missing
- pandas support
- Adding new algorithms

If you are really interested in going down this path, please read on:


-----------------
Project Structure
-----------------

Each fitting algorithm resides in its own module. All the functions get
imported into the base :mod:`skg` namespace. Each module should contain a
function called `model` that applies the fitting parameters to a given set of
x-values, either raveled or along a particular axis (assuming the function is
1D). Multiple algorithms that fit to the same model can live in the same
module.


-------
Testing
-------

Testing is done using the :ref:`ptyest` framework. A test module for every main
module exists in the :mod:`tests` package.

Tests for new modules are generated in a semi-automated manner (still WIP). All
the modules containing a fitting function and a model will be tested against
randomly generated inputs, and checked for speed and quality. The quality of
each algorithm will be assessed based on these tests. Quality has three
categories: speed, accuracy and usefulness.

- Speed is a benchmark against :func:`scipy.optimize.curve_fit`. An algorithm
  that is slower than a non-linear optimizer starting with default parameters
  is not deemed very useful.
- Accuracy is checked by making sure that the fit is within reasonable bounds
  of the values computed by :func:`scipy.optimize.curve_fit`. Reasonableness is
  a function of the analytically derived partial derivatives of the model with
  respect to the parameters.
- Usefulness is a measure of how many iterations
  :func:`scipy.optimize.curve_fit` saves by using the algorithm as an initial
  guess. Another informal metric is the combined runtime of the algorithm and
  :func:`~scipy.optimize.curve_fit` vs. the runtime of just
  :func:`~scipy.optimize.curve_fit` with default parameters. If the latter
  exceeds the former, that's a win.
