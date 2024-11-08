# -*- encoding: utf-8 -*-
import os
import toml

pyproject_path = os.path.join(os.path.join('..','pyproject.toml'))
with open(pyproject_path, 'r') as f:
    pyproject_data = toml.load(f)

__version__ = pyproject_data['tool']['poetry']['version']
__author__ = pyproject_data['tool']['poetry']['authors']
__original_url__ = "https://github.com/zonca/planck-levelS/tree/master/Beam"
