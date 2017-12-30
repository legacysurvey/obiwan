=======
Obiwan
=======

Obiwan is a Monte Carlo method for adding fake galaxies and stars to the Legacy Survey imaging data. We forward model astrophysical sources in the images with the legacypipe code, by identifying Signal to Noise (S/N) > 6 sources and finding the best fitting 2D model for a star or galaxy by minimizing a regularized L2 Loss function.

Use the links at the top of the page to learn more about the code or see the table of contents below.

.. toctree::
   :caption: Table of Contents

Using Obiwan at the National Energy Researc Scientific Computing Center (NERSC)

* :doc:`Production Runs <howto/PRODUCTION_RUN>`
* :doc:`Data Model <howto/OUTPUTS>`

The following science is being done

:ref:`Travis continuous integration tests of the entire pipeline <nb/TestCases.ipynb>`
----------------------------------------------------------------------------------------

The module :mod:`tests.end_to_end.test_datasets` carries out these tests. See the visualizations below for what these tests do. 

* :ref:`Stars <nb/TestCases.ipynb#star:-testcase_DR_z>`
* :ref:`Galaxies <nb/TestCases.ipynb#allblobs=False:-testcase_DR_z>`
* :ref:`Sources on the edges of images <nb/TestCases.ipynb#onedge:-testcase_DR_z>`
* :ref:`Early coadds <nb/TestCases.ipynb#Early-coadds:-elg-grz>`

Changelog
----------

* :doc:`changes`

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


