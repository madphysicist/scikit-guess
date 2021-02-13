===============
Getting Started
===============


------------
Installation
------------

Installation through pypi is preferred, e.g. with ::

    $ pip install scikit-guess

Installation from source is also an option. The latest version can be obtained
from GitHub via ::

    $ git clone github.com:madphysicist/scikit-guess

You can then install in a couple of different ways. The standard is ::

    $ python setup.py install

To install for development purposes, or just to have the latest bleeding-edge
version always running ::

    $ python setup.py develop


------------
Dependencies
------------

This project strives to be a true scikit, and limit its runtime dependencies to
just `numpy`_ and `scipy`_. `pandas`_ support will be added at some point, but
may not result in a dependency, certainly not a mandatory one. Numpy must be
version 1.7 or later, since that is when the ``where`` keyword was introduced
to ufuncs.

At this stage, the code is written in pure python, with all the extensions
being implemented through the dependencies. That may change at some point in
the future.

Python 3.6+ is supported.

Testing is done with `pytest`_ and has an optional dependency on `matplotlib`_
for debugging images. Another optional dependency is `pytest-pep8`_.

Building the documentation requires `sphinx`_ (version >= 1.8 preferred, see
:ref:`below <start-docs>`). The  `Sphinx RTD Theme`_ is an optional
dependency. The theme will fail over to the built-in `Alabaster Theme`_ instead
if RTD is missing.


-----
Tests
-----

To validate your install, you can run the tests. To test an installed version::

    >>> import skg
    >>> skg.test()

To test a source distribution, you can use either ::

    $ python setup.py test

or ::

    $ pytest

The latter is preferred.

All the usual `pytest` command line options are allowed. To enable plotting of
some of the fixtures and output results, a ``--plots`` flag is provided::

    >>> skg.test('--plots')

or ::

    $ python setup.py test --pytest-args=--plots

or ::

    $ pytest --plots


.. _start-docs:

-------------
Documentation
-------------

To build the documentation from a source install, run either ::

    $ python setup.py doc

or ::

    $ make -C doc clean html

The former should work on all systems. The latter will probably not work on
Windows. Instead, you will have to switch to the `doc` directory and run
without the ``-C doc`` option.

Currently, both versions run with nitpicky mode by default. The `setup.py`
version requires `sphinx`_ >= 1.8 for this to work. If you would like to use an
older version, comment out the ``nitpicky`` option in `setup.cfg`.

The documentation output goes to `build/doc` in the root directory with both
versions.

The HTML documentation uses a small amount of CSS3, and will likely not render
100% accurately on IE8 or older. The discrepancies should be quite minimal
however.

It is also possible to build the documentation as LaTeX using ::

    $ make -C doc clean latex

Direct output to PDF through the ``pdf`` or ``pdflatex`` targets is not
supported. However, the LaTeX build can be used to genrate PDFs::

    $ make -C build/doc/latex pdf-all

This requires some additional dependencies, like a LaTeX system and `latexmk`.


.. include:: /link-defs.rst
