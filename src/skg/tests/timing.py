"""
This file is not part of the formal functional tests for the package.
Instead, it contains a set of functions that time various ways of
computing the datasets required by the scikit. The timing is computed
using the :mod:`timeit` module.

Each function times the methods as fairly as possible against each
other, and recommends the one to use. The timing on the author's
original development machine dictated the layout of many of the methods
in the scikit.

This file can be run as a script on any given target machine to optimize
for a particular architecture. Modifications to the source code can be
made accordingly. The usual regression tests should work regardless of
the method selected.
"""

import timeit


def run_all():
    pass


if __name__ == '__main__':
    run_all()
