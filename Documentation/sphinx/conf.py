# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import sys, os
from pathlib import Path

sys.path.insert(0, str(Path('..', '..').resolve()))
sys.path.insert(0, str(Path('..', '..', 'scripts').resolve()))
print(str(Path('..','..').resolve()))
print(str(Path('..','..', 'scripts').resolve()))
project = 'analysator'
copyright = '2024, University of Helsinki'
author = 'Sameli'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc'
    ]

templates_path = ['_templates']
# Create a dummy _templates folder if it does not exist
try:
    os.mkdir("_templates")
except:
    pass


exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_logo = "logo_color.png"
html_static_path = ['_static']

# Create a dummy _static folder if it does not exist
try:
    os.mkdir("_static")
except:
    pass
