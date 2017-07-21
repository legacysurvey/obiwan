# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
============
obiwan
============

A package for injecting artificial stars and galaxy-targets (ELGs, LRGs, and QSOs) into the individual images of the DESI Legacy Surveys, DECALS, MzLS, and BASS. The legacypipe/Tractor pipeline is rerun on the modified image data, and the resulting *simulated* Tractor catalogues form a Data Release (DR) resembling a DR of real objects. 
"""
#
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
# The line above will help with 2to3 support.
#
# Set version string.
#
from ._version import __version__
