|
|

=======
Obiwan
=======

|
|

.. title:: obiwan docs

.. toctree::
   :caption: Table of Contents

A Monte Carlo method for adding fake galaxies and stars to images from the Legacy Survey and re-processing the modified images with our pipeline. Our pipeline forward models galaxies and stars in the multi-color images by detecting sources with Signal to Noise (S/N) greater than 6 and minimizing the regularized L2 Loss function for various models for the shapes of stars and galaxies.


Why the name "obiwan"?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Just as Obi-Wan Kenobi was the "only hope" in Star Wars: Episode IV - A New Hope; **obiwan**, by virtue of its Monte Carlo method, is possibly the only hope for removing all systematics in the sample of galaxies we select from the imaging data with our `Legacypipe pipeline <https://github.com/legacysurvey/legacypipe>`_. These are commonly referred to as "imaging systematics" since they are related to image quality, the telescope that took the image, and the bias, bugs, completeness, etc. of the `Legacypipe pipeline <https://github.com/legacysurvey/legacypipe>`_ itself. The Sloan Digital Sky Survey (SDSS), showed that these "imaging systematics" could be removed by measuring correlations between the number of galaxies and image quality (namely, stellar density, seeing, galactic extinction, sky background, and photometric offsets). The situation is far worse for the Legacy Surveys because our images come from three telescopes (SDSS used one), and we take exposures of the same part of the sky, in multiple bands, over timescales of years (SDSS imaged each part of the sky in all bands nearly simultaneously, e.g. over a few minutes). In addition, the cosmological signal we are looking for in the images, about **FILL IN** arcsec, is the size of the field of view of our telescopes' cameras (the SDSS telescope had a field of view of **FILL IN**).

It's probably incorrect to describe our Monte Carlo additions of fake galaxies to the imaging data as the "only hope" of identifying imaging systematics, but it should be one of the most robust and self consistent methods for doing so.

End-to-End test cases
^^^^^^^^^^^^^^^^^^^^^^

We created various tests that add fake galaxies to small 200x200 pixel multi-color images, run our pipeline, and assert that the fake galaxies are detected and have shape and multi-color brightness very close to truth.

The Turorials section shows what these End-to-End tests look like, which come from this :ref:`ipython notebook <nb/TestCases.ipynb>`. We carry out the tests with :mod:`tests.end_to_end.test_datasets` python module.

Running the code at NERSC
^^^^^^^^^^^^^^^^^^^^^^^^^^

Using Obiwan at the National Energy Researc Scientific Computing Center (NERSC)

* :doc:`Production Runs <howto/PRODUCTION_RUN>`
* :doc:`Data Model <howto/OUTPUTS>`

Acknowledgements
^^^^^^^^^^^^^^^^^^

See the `offical Acknowledgements <http://legacysurvey.org/#Acknowledgements>`_ for the Legacy Survey.


Changelog
^^^^^^^^^^^

* :doc:`changes`

Indices and tables
^^^^^^^^^^^^^^^^^^^

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


