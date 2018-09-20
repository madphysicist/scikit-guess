===============
Getting Started
===============


------------
Installation
------------

Installation through pypi is preferred, e.g. with ::

    pip install scikit-guess

Installation from source is also an option. The latest version can be obtained
from GitHub via ::

    git clone github.com:madphysicist/scikit-guess

You can then install in a couple of different ways. The standard is ::

    python setup.py install

To install for development purposes, or just to have the latest bleeding-edge
version always running ::

    python setup.py develop


------------
Dependencies
------------

This project strives to be a true scikit, and limit its runtime dependencies to
just `numpy`_ and `scipy`_. `pandas`_ support will be added at some point, but
may not result in a dependency, certainly not a mandatory one.

At this stage, the code is written in pure python, with all the extensions
being implemented through the dependencies. That may change at some point in
the future.

Python 2.7 and 3.4+ are supported.

Testing is done with `pytest`_ and has an optional dependency on `matplotlib`_
for debugging images. Another optional dependency is `pytest-pep8`_.

Building the documentation requires `sphinx`_. The  `Sphinx RTD Theme`_ is an
optional dependency. The theme will fail over to the built-in
`Alabaster Theme`_ instead if RTD is missing.


-------
Testing
-------

To validate your install, you can run the tests. To test an installed version::

    >>> import skg
    >>> skg.test()

To test a source distribution, you can either use ::

    $ python setup.py test

or ::

    $ pytest

The latter is preferred.


.. include:: /link-defs.rst