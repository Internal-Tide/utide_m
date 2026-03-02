import os
import sys
sys.path.insert(0, os.path.abspath('..'))

project = 'UTide'
copyright = '2026, UTide Contributors'
author = 'UTide Contributors'
release = '1.0'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.mathjax',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

language = 'zh_CN'

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

autodoc_member_order = 'bysource'
napoleon_google_docstring = False
napoleon_use_param = True
napoleon_use_ivar = True
