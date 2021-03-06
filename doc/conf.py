# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
import os
import sys
sys.path.insert(0, os.path.abspath('../src'))
import skg.version


# -- Project information -----------------------------------------------------

project = 'scikit-guess'
copyright = '2021, Joseph R. Fox-Rabinovitz'
author = 'Joseph R. Fox-Rabinovitz'

# The short X.Y version
version = ''
# The full version, including alpha/beta/rc tags
release = skg.version.__version__


# -- General configuration ---------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.coverage',
    'sphinx.ext.doctest',
    'sphinx.ext.ifconfig',
    'sphinx.ext.imgmath',
    'sphinx.ext.intersphinx',
    'sphinx.ext.napoleon',
    'sphinx.ext.todo',
    'sphinx.ext.viewcode',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = None

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path .
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# Don't add () after function names
add_function_parentheses = False

# Number figures and tables
numfig = True


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
try:
    import sphinx_rtd_theme
except ImportError:
    html_theme = 'alabaster'
else:
    html_theme = 'sphinx_rtd_theme'

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
# html_theme_options = {}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Custom sidebar templates, must be a dictionary that maps document names
# to template names.
#
# The default sidebars (for documents that don't match any pattern) are
# defined by theme itself.  Builtin themes are using these templates by
# default: ``['localtoc.html', 'relations.html', 'sourcelink.html',
# 'searchbox.html']``.
#
# html_sidebars = {}


# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'scikit-guessdoc'


# -- Options for LaTeX output ------------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',

    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, 'scikit-guess.tex', 'scikit-guess Documentation',
     'Joseph R. Fox-Rabinovitz', 'manual'),
]


# -- Options for manual page output ------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'scikit-guess', 'scikit-guess Documentation',
     [author], 1)
]


# -- Options for Texinfo output ----------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, 'scikit-guess', 'scikit-guess Documentation',
     author, 'scikit-guess', 'One line description of project.',
     'Miscellaneous'),
]


# -- Extension configuration -------------------------------------------------

# -- Options for autosummary extension ---------------------------------------

# Auto-generate stubs
autosummary_generate = True

# -- Options for intersphinx extension ---------------------------------------

# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {
    'lmfit': ('https://lmfit.github.io/lmfit-py/', None),
    'matplotlib': ('https://matplotlib.org/', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'pytest': ('https://docs.pytest.org/en/latest/', None),
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/reference/', None),
}

# -- Options for napoleon extension ------------------------------------------

napoleon_use_param = True
napoleon_type_aliases = {
    'array-like': ':term:`array-like <numpy:array_like>`',
    'array_like': ':term:`numpy:array_like`',
}

# -- Options for todo extension ----------------------------------------------

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = True


# -- Project configuration ---------------------------------------------------

# Content generators
content_preprocess = ['appendices/reei']

def preprocess_content():
    """
    This is a hack to ensure that the appropriate content is generated
    for the figures and tables in the paper. It checks the folders
    listed in `content_preprocess` and runs each file there.

    This is not a really good stand-in for say a proper directive to
    auto-generate the content right there and then into a configurable
    folder, etc., but it will have to do.

    If such an extension were to be created, it would add directives for
    `auto-figure` and `auto-table` as a minimum. These would extend the
    normal `figure` and `table` directives to avoid hassle. They would
    just add a mixin that would add a couple of additional configuration
    arguments:

    1. `cache`: "yes" or "no" value. If "yes", check the output file for
       date before running the generator. If the file exists and is
       newer than the script file, use it as-is. If "no", always run the
       script.
    2. `generator`: The name of a python module, followed by `::` and a
       function name. The default function name is generate_figure or
       generate_table, as the case may be. The whole module is imported,
       and the selected function is run with the name of the output file
       as an input.
    3. `source`: `auto-figure` would already require a file name to be
       supplied, but `auto-table` would not. `auto-table` would
       therefore require a `source` option to supply that. The default
       could be to just replace itself with the text of the table and
       go for another round of parsing instead of including a file.
    """
    from glob import glob
    from os import getcwd
    from os.path import abspath, basename, join, splitext
    from contextlib import contextmanager

    @contextmanager
    def path_context(location):
        sys.path.insert(0, abspath(join(getcwd(), location)))
        yield
        del sys.path[0]

    with path_context('..'):
        from setup import import_file

    for folder in globals().get('content_preprocess', []):
        for file in glob(join(folder, '*.py')):
            name = splitext(basename(file))[0]
            import_file(name, file)


def setup(app):
    """
    Add custom stylesheet(s).
    """
    app.add_css_file('css/custom.css')
    preprocess_content()
