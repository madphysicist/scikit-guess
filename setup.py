#!/usr/bin/env python

from __future__ import absolute_import

from os.path import dirname, join
import sys

from setuptools import setup
from setuptools.command.test import test as TestCommand


DIST_NAME = 'scikit-guess'

LICENSE = 'BSD 2-Clause License'
DESCRIPTION = 'Non-iterative initial parameter guesses for fitting routines'

AUTHOR = 'Joseph R. Fox-Rabinovitz'
AUTHOR_EMAIL = 'joseph.r.fox-rabinovitz@nasa.gov'

MAINTAINER = 'Joseph R. Fox-Rabinovitz'
MAINTAINER_EMAIL = 'jfoxrabinovitz@gmail.com'

CLASSIFIERS = (
    'Development Status :: 2 - Pre-Alpha',
    'Intended Audience :: Developers',
    'Intended Audience :: Education',
    'Intended Audience :: Financial and Insurance Industry',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: BSD License',
    #'Operating System :: Microsoft :: Windows',
    'Operating System :: POSIX :: Linux',
    #'Operating System :: Unix',
    #'Operating System :: MacOS',
    'Programming Language :: Python',
    'Programming Language :: Python :: 2',
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: Implementation :: CPython',
    #'Programming Language :: Python :: Implementation :: PyPy',
    #'Programming Language :: Python :: Implementation :: Jython',
    #'Programming Language :: Python :: Implementation :: IronPython',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Information Analysis',
    'Topic :: Scientific/Engineering :: Mathematics',
    'Topic :: Office/Business :: Financial :: Investment',
)


COMMANDS = {}


class PyTest(TestCommand):
    """
    Suggested by pytest documentation to avoid dependency on
    `pytest-runner` when integrating with `test` command.
    """
    user_options = [('pytest-args=', 'a', 'Arguments to pass to pytest')]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = ''

    def run_tests(self):
        import shlex

        # import here, cause outside the eggs aren't loaded
        import pytest

        args = shlex.split(self.pytest_args)
        args.insert(0, '-p')
        args.insert(1, 'skg.tests.options')
        errno = pytest.main(args)
        sys.exit(errno)

COMMANDS['test'] = PyTest


try:
    from sphinx.setup_command import BuildDoc
    COMMANDS['build_sphinx'] = BuildDoc
except ImportError:
    pass


def import_file(name, location):
    """
    Imports the specified python file as a module, without explicitly
    registering it to `sys.modules`.
    """
    if sys.version_info[0] == 2:
        # Python 2.7-
        from imp import load_source
        mod = load_source(name, location)
    elif sys.version_info < (3, 5, 0):
        # Python 3.4-
        from importlib.machinery import SourceFileLoader
        mod = SourceFileLoader(name, location).load_module()
    else:
        # Python 3.5+
        from importlib.util import spec_from_file_location, module_from_spec
        spec = spec_from_file_location(name, location)
        mod = module_from_spec(spec)
        spec.loader.exec_module(mod)
    return mod


def version_info():
    """
    Jump through some hoops to import version.py for the different
    versions of Python.

    https://stackoverflow.com/a/67692/2988730
    """
    location = join(dirname(__file__) or '.', 'src', 'skg', 'version.py')
    mod = import_file('version', location)
    return mod.__version__


def long_description():
    """
    Reads in the README.md and CHANGELOG.md files, separated by two
    newlines.
    """
    with open('README.md') as readme, open('CHANGELOG.md') as changes:
        return '%s\n\n%s' % (readme.read(), changes.read())


if __name__ == '__main__':
    setup(
        name=DIST_NAME,
        version=version_info(),
        license=LICENSE,
        description=DESCRIPTION,
        long_description_content_type='text/markdown',
        long_description=long_description(),
        author=AUTHOR,
        author_email=AUTHOR_EMAIL,
        maintainer=MAINTAINER,
        maintainer_email=MAINTAINER_EMAIL,
        classifiers=CLASSIFIERS,
        url='https://github.com/madphysicist/scikit-guess',
        project_urls={
            'Bugs': 'https://github.com/madphysicist/scikit-guess/issues',
            'Documentation': 'https://scikit-guess.readthedocs.io/en/latest/',
        },
        packages=['skg', 'skg.tests'],
        package_dir={'': 'src'},
        install_requires=[
            'numpy >= 1.7',
            'scipy',
        ],
        tests_require=['pytest'],
        cmdclass=COMMANDS,
        extras_require={
            'pandas': ['pandas'],
            'test-plots': ['matplotlib'],
            'test-pep8': ['pytest-pep8'],
            # TODO: Some of the sphinx extensions may need to go in here.
            'docs': [
                'matplotlib',
                'numpy',
                'scipy',
                'sphinx >= 1.8'
            ],
            'docs-rtd': ['sphinx_rtd_theme'],
        }
    )
