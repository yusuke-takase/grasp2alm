# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
import toml

module_path = os.path.abspath(os.path.join('..','..','grasp2alm'))
pyproject_path = os.path.abspath(os.path.join('..','..','pyproject.toml'))
sys.path.insert(0, module_path)

with open(pyproject_path, 'r') as f:
    pyproject_data = toml.load(f)

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'grasp2alm'
copyright = '2024, Yusuke Takase'
author = pyproject_data['tool']['poetry']['authors']
release = pyproject_data['tool']['poetry']['version']

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
    'nbsphinx',
    'sphinx.ext.autosectionlabel',
    'pydata_sphinx_theme'
]

autosummary_generate = True
autosummary_imported_members = True
autosectionlabel_prefix_document = True
autoclass_content = "class"

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

#html_theme = 'alabaster'
html_theme = 'pydata_sphinx_theme'
#html_static_path = ['_static']
