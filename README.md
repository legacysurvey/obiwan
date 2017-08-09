# obiwan
An extension of the legacypipe/Tractor pipeline for the [DECALS, MzLS, and BASS Legacy Surveys](http://legacysurvey.org/dr4/description). This package injects artificial stars and galaxy-targets (ELGs, LRGs, and QSOs) into individual images for uncerstanding imaging systematic for [DESI](https://desi.lbl.gov). It re-runs the legacypipe/Tractor pipeline on the modified images and the resulting _simulated_ Tractor catalogues can be used as a Data Release.

[![Docs](https://readthedocs.org/projects/obiwan/badge/?version=latest)](http://obiwan.readthedocs.org/en/latest/)
[![Build Status](https://travis-ci.org/legacysurvey/obiwan.png)](https://travis-ci.org/legacysurvey/obiwan)
[![Coverage Status](https://coveralls.io/repos/github/legacysurvey/obiwan/badge.svg?branch=master)](https://coveralls.io/github/legacysurvey/obiwan)

# Documentation
------------------

Please visit [obiwan on locally built Docs](https://legacysurvey.github.io/obiwan)

The Docs are built locally from these [instructions](https://github.com/legacysurvey/obiwan/tree/gh-pages/README.md), because the galsim install is so complicated that getting .travis.yml installing everything correctly I didn't want to repeat with the Sphinx/Read the Docs build script. The deprecated documentation is here [obiwan on Read the Docs](http://obiwan.readthedocs.org/en/latest/)

License
=======

obiwan is free software licensed under a 3-clause BSD-style license. For details see
the ``LICENSE.rst`` file.
