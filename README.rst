===========
obiwan
=========

An extension of the legacypipe/Tractor pipeline for the `DECALS, MzLS, and BASS Legacy Surveys <http://legacysurvey.org/dr4/description>`_. This package injects artificial stars and galaxy-targets (ELGs, LRGs, and QSOs) into individual images for uncerstanding imaging systematic for `DESI <https://desi.lbl.gov>`. It re-runs the legacypipe/Tractor pipeline on the modified images and the resulting *simulated* Tractor catalogues can be used as a Data Release.

|docs| |build-status| |coverage|

.. |docs| image:: https://readthedocs.org/projects/obiwan/badge/?version=latest
    :alt: Documentation Status
    :scale: 100%
    :target: http://obiwan.readthedocs.io/en/latest/?badge=latest

.. |build-status| image:: https://travis-ci.org/legacysurvey/obiwan.svg?branch=master
    :alt: Build Status
    :scale: 100%
    :target: https://travis-ci.org/legacysurvey/obiwan

.. |coverage| image:: https://coveralls.io/repos/github/legacysurvey/obiwan/badge.svg?branch=master
    :alt: Coverage Status
    :scale: 100%
    :target: https://coveralls.io/github/legacysurvey/obiwan


Documentation
^^^^^^^^^^^^^^

See our `RTD page <http://obiwan.readthedocs.io/en/latest/?badge=latest>

Run at `NERSC <http://www.nersc.gov/>`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

See the :doc:`Run Instructions <bin/run_atnersc/README.md>`.

Need Help?
^^^^^^^^^^^^^

Email us: desi-image-sims 'at' googlegroups.com 

Acknowledgements
^^^^^^^^^^^^^^^^

See the `Legacy Survey Acknowledgements` <http://legacysurvey.org/#Acknowledgements>`.

License
^^^^^^^^^^^

obiwan is free software licensed under a 3-clause BSD-style license. For details see
our :doc:`LICENSE <LICENSE.rst>`.
