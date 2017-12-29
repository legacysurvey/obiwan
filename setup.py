#!/usr/bin/env python
# -*- encoding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function

import io
import re
from glob import glob
import os

from setuptools import find_packages
from setuptools import setup


setup_kwargs= dict(
    name='obiwan',
    version='1.2',
    license='BSD',
    description='Code to Monte Carlo fake sources through the Legacypipe pipeline',
    long_description='',
    author='Kaylan Burleigh, John Moustakas',
    author_email='kburleigh@lbl.gov',
    url='https://github.com/legacysurvey/obiwan',
    packages=find_packages('py/obiwan'),
    package_dir={'': 'py/obiwan'},
    py_modules=[os.path.splitext(os.path.basename(path))[0] 
                for path in glob('py/obiwan/*.py')],
)
with open('README.md') as readme:
    setup_kwargs['long_description'] = readme.read()

setup(**setup_kwargs)



