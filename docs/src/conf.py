# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = "i-PI"
copyright = "2021, The i-PI developers"
author = "The i-PI developers"

# The full version, including alpha/beta/rc tags
import configparser

config = configparser.ConfigParser()
config.read("../../setup.cfg")
release = config["metadata"]["version"]

# -- General configuration ---------------------------------------------------
needs_sphinx = "3.2"

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ["sphinxcontrib.bibtex", "sphinx.ext.todo"]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**/input_ref_sections/*"]


bibtex_bibfiles = ["references.bib"]
# bibtex_encoding = 'latin'
# bibtex_default_style = 'unsrt'

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_title = "i-PI documentation pages"
html_theme = "furo"
html_theme_options = {
    "top_of_page_buttons": [],
}
html_favicon = "../_static/favicon-ipi.png"
html_logo = "../figures/ipi-logo.svg"

html_css_files = [
    "custom_styles.css",
]

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["../_static"]

# -- Options for LaTeX output -------------------------------------------------

latex_theme = "howto"

latex_elements = {"fontenc": "\\usepackage[LGR,T1]{fontenc}"}
