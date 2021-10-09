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
import os
import io
import sys
import shutil
# sys.path.insert(0, os.path.abspath('.'))




# -- Project information -----------------------------------------------------

project = 'Demuxafy'
copyright = '2021, Drew Neavin'
author = 'Drew Neavin'

# The full version, including alpha/beta/rc tags
release = '0.0.1'

# XXX: POSSIBLY KIND OF A HACK
# import httpolice
# import httpolice.inputs
# import httpolice.reports.html
 
# Use hack developed by https://github.com/vfaronov/httpolice/blob/5bf58456b37d2f073094e6f39f168b41fda5ef3d/doc/conf.py#L33-L63
# to import html generated for calculator as it's own page
# ``test.html`` with an example HTML report.

 
# But for that to work, we need to actually build those pages first, and put
# them into the ``html_extra_path``. I can't think of a more appropriate place
# to trigger their build than right here. After all, Sphinx's docs do say that
# ``conf.py`` "can execute arbitrarily complex code", so maybe it's OK.

# See also: https://stackoverflow.com/q/38547509/200445

# if os.path.exists('_extra'):
#     shutil.rmtree('_extra')
# os.mkdir('_extra')

# with io.open('_extra/test.html', 'wb') as notices_file:
#     httpolice.reports.html.list_notices(notices_file)

# with io.open('_extra/showcase.html', 'wb') as showcase_file:
#     exchanges = list(httpolice.inputs.combined_input(
#         ['../test/combined_data/showcase.https']))
#     for exch in exchanges:
#         httpolice.check_exchange(exch)
#     httpolice.html_report(exchanges, showcase_file)

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc',
				'sphinx_tabs.tabs',
				'sphinx_copybutton',
				'sphinxemoji.sphinxemoji',
				'sphinx_togglebutton']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']


# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = [".rst"]

# The master toctree document.
master_doc = "index"

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "colorful"



# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# html_theme = 'alabaster'
html_theme = "sphinx_rtd_theme"


# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
html_theme_options = {
    "navigation_depth": 1,
    # "logo_only": True,
}

# html_logo = "logo.png"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_extra_path = ['_extra']

html_css_files = [
    'customs.css',
] 