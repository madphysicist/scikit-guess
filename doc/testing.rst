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

    conda create --name skg-testing-py3.6 --no-default-packages python=3.6 nomkl numpy scipy pytest sphinx sphinx_rtd_theme
    conda create --name skg-testing-py3.7 --no-default-packages python=3.7 nomkl numpy scipy pytest sphinx sphinx_rtd_theme
    ...
