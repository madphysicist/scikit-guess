#!/bin/bash

# This is a small script that encapsulates commands for building and
# uploading the distros. To run, create virtual environments for any
# python 3 interpreter under ../venvs/skg and make sure that the
# wheel package is installed for both.

source ../.venv/skg/bin/activate
python setup.py sdist bdist_wheel
make -C doc clean html
twine upload --repository-url https://test.pypi.org/legacy/ dist/*
deactivate
