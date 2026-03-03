==========================
Building documentation
==========================

Documentation is automatically built and deployed using Sphinx as part of the continuous integration of the repository. If you want to build the documentation locally, you need to install a few more packages::

   python -m venv CI_env
   . CI_env/bin/activate
   pip install ../analysator[all]
   pip install sphinx numba yt geopack sphinx-rtd-theme

after which you can generate e.g. the html documentation::

   cd Documentation/sphinx
   make dirhtml

