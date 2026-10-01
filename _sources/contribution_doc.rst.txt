Contributing to documentation
=============================

Guiding principles
------------------

We recommend documentation contributors to familiarise themselves with the [Diátaxis](https://diataxis.fr/) framework when deciding which type of documentation to write. The key approach is to distinguish between four fundamentally different categories of documentation: tutorials, how-to guides, explanation material, and reference material. But most importantly, the idea is to work on documentation piece by piece whenever the occasion presents itself. If you find a mistake in the documentation, if you write new code and its documentation, if you contribute a tutorial – all of this is welcome in any stage of readiness! **A contribution to improve the documentation, however small, is better than no contribution at all.**


Automatic building and deployment
---------------------------------

This documentation is generated via Sphinx-autodoc from the Analysator source. This page is a work in progress, with docstrings and the documentation layout in flux. Contributions via the Analysator repository are welcome.

Add your content using ``reStructuredText`` syntax. See the
`reStructuredText <https://www.sphinx-doc.org/en/master/usage/restructuredtext/index.html>`_
documentation for details.

Documentation is automatically built and deployed using Sphinx as part of the continuous integration of the repository.


Building documentation
----------------------

When you create new documentation, before posting a pull request, please build the documentation locally and check the results.

In order to build the documentation, you need to install a few packages::

   python -m venv CI_env
   . CI_env/bin/activate
   pip install ../analysator[all]
   pip install sphinx numba yt geopack sphinx-rtd-theme

after which you can generate e.g. the html documentation::

   cd Documentation/sphinx
   make dirhtml

and open it in your favourite browser.


Sphinx coverage stats
---------------------

.. include:: coverage.rst
   :literal:

