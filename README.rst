===========
obiwan
===========

**Obiwan** is a Monte Carlo method for adding fake galaxies and stars to images from the Legacy Survey and re-processing the modified images with our `Legacysurvey/Tractor pipeline <https://github.com/legacysurvey/legacypipe>`_. The pipeline forward models galaxies and stars in the multi-color images by detecting sources with Signal to Noise (S/N) greater than 6 and minimizing the regularized L2 Loss function for various models for the shapes of stars and galaxies.

Build Status
^^^^^^^^^^^^^

|build-status| |coverage|

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

Our documentation is hosted on GitHub Pages, `Obiwan Docs <https://legacysurvey.github.io/obiwan/>`_.

Contact Us
^^^^^^^^^^^

* Email: desi-image-sims 'at' googlegroups.com 

^^^^^^^^^^^^^^^^^^

See the `offical acknowledgements <http://legacysurvey.org/#Acknowledgements>`_ for the Legacy Survey.

License
^^^^^^^^^^^

obiwan is free software licensed under a 3-clause BSD-style license. For details see the `LICENSE <https://github.com/legacysurvey/obiwan/blob/master/LICENSE.rst>`_.
