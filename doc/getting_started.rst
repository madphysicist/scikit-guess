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

This project strives to be a true scikit, and limit it's dependencies to just
numpy and scipy. Pandas support will be added at some point, but may not result
in a dependency, certainly not a mandatory one.

At this stage, the code is written in pure python, with all the extensions
being implemented through the dependencies. That may change at some point in
the future.

Python 2.7 and 3.4+ are supported.

Building the documentation requires :ref:`sphinx`.
