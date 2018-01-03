===========
obiwan
===========

**Obiwan** is a Monte Carlo method for adding fake galaxies and stars to images from the Legacy Survey and re-processing the modified images with our `Legacysurvey/Tractor pipeline <https://github.com/legacysurvey/legacypipe>`_. The pipeline forward models galaxies and stars in the multi-color images by detecting sources with Signal to Noise (S/N) greater than 6 and minimizing the regularized L2 Loss function for various models for the shapes of stars and galaxies.

Build Status
^^^^^^^^^^^^^

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

Our documentation is hosted on `ReadTheDocs <http://obiwan.readthedocs.io/en/latest/?badge=latest>`_.

How to run the code
^^^^^^^^^^^^^^^^^^^^

Detailed instructions for running obiwan on the National Energy Research Scientific Computing Center (NERSC) supercomputers are below.

* `How to Run at Data Release <https://github.com/legacysurvey/obiwan/blob/master/doc/howto/PRODUCTION_RUN.rst>`_
* `Description of the outputs <https://github.com/legacysurvey/obiwan/blob/master/doc/howto/OUTPUTS.rst>`_
* `How to Train a CNN <https://github.com/legacysurvey/obiwan/blob/master/doc/deeplearning.rst>`_

Contact Us
^^^^^^^^^^^

* Email: desi-image-sims 'at' googlegroups.com 

^^^^^^^^^^^^^^^^^^

See the `offical acknowledgements <http://legacysurvey.org/#Acknowledgements>`_ for the Legacy Survey.

License
^^^^^^^^^^^

obiwan is free software licensed under a 3-clause BSD-style license. For details see the `LICENSE <https://github.com/legacysurvey/obiwan/blob/master/LICENSE.rst>`_.