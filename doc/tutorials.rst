.. _tutorials:

Tutorials
==========

Each section below shows lessons learned via ipython notebooks and/or exploratory data analysis.

End-to-End test cases
-----------------------

We created our End-to-End test cases using the `create_testcase.py <https://github.com/legacysurvey/legacypipe/blob/master/py/legacyanalysis/create_testcase.py>`_ script in our Legacypipe pipeline. You can visualize the various Test Cases by clicking on the links below.

.. end-to-end-tests

.. toctree::
  :maxdepth: 3

  nb/TestCases.ipynb

Modeling the color, shape, and redshift of (ELG-like) galaxies
--------------------------------------------------------------------

.. toctree::
  :maxdepth: 3

  nb/eBOSS_DESI_Priors.ipynb

Deep Learning 
--------------
Obiwan creates millions of images of fake galaxies. If they are representative fo real galaxies, a Neural Network should not be able to tell them apart. Read about the project :doc:`here <deeplearning>`. 

Example code is below.

.. toctree::
  :maxdepth: 3

  nb/CNN.ipynb


Test Region: 10 deg2
---------------------

The first true test of Obiwan was a 10 deg2 test region, described here :ref:`PRODUCTION_RUN.md`. 

.. toctree::
  :maxdepth: 3

  nb/10deg2_region.ipynb


Correlation Functions
---------------------

We use `nbodykit` to compute correlation functions. The randoms are Obiwan Randoms.

.. toctree::
  :maxdepth: 3

  nb/AngCorr_Func.ipynb

.. _nbodykit: http://github.com/bccp/nbodykit









