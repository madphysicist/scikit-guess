Testing SKG
===========

.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. automodule:: skg.tests


----------------
Additional Notes
----------------

Test environments can be created under :program:`conda` with the following
commands::

    conda create --name skg-testing-py2.7 --no-default-packages python=2.7 nomkl numpy scipy pytest sphinx sphinx_rtd_theme
    conda create --name skg-testing-py3.4 --no-default-packages python=3.4 nomkl numpy scipy pytest sphinx sphinx_rtd_theme
    ...
