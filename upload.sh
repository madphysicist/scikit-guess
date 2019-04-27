#!/bin/bash

# This is a small script that encapsulates commands for building and
# uploading the distros. To run, create virtual environments for any
# python 2 and 3 interpreter under ./venvs and make sure that the
# wheel package is installed for both.

source ./venvs/2.7/bin/activate
python setup.py sdist bdist_wheel
source ./venvs/3.6/bin/activate
python setup.py sdist bdist_wheel
deactivate
make -C doc clean html

twine upload --repository-url https://test.pypi.org/legacy/ dist/*
