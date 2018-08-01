******************
Post Processing
******************

All post processing of an obiwan production run is done by ``obiwan/py/obiwan/runmanager/derived_tables.py``
A single fits table is created per brick, which I call a "derived table". It contains the randoms table ra, dec, fluxes, and shapes, fluxes and shapes actually added to the images, and the tractor measurements (if detected) for each of these. A few bit masks are created, one says which injected sources were recovered and modeled by legacypipe, which of those are thought to be coincident with real galaxies from DR3 or DR5 etc. Another bit mask says which of the injected sources would pass target selection based on their true fluxes and which pass based on their tractor measured fluxes.

Takes a list of bricks and creates each table in an embarrassingly parallel fashion using mpi4py.

There are two modes: ``randoms`` and ``summary``, randoms is the derived table while summary is a table containing various stats about each brick, e.g., number of sources injected, average depth of sources, fraction of injected sources detected by legacypipe.

Run it as a batch job using this script slurm_derived_tables.sh_

.. _slurm_derived_tables.sh: https://github.com/legacysurvey/obiwan/blob/master/bin/slurm_derived_tables.sh

The per-brick tables can be combined into a single table using  ``obiwan/py/obiwan/runmanager/merged_tables.py``. There are two modes: ``parallel`` and ``serial``. Parallel is run first and it combines the per-brick tables into < 100 tables (a much easier number than > 10,000). Serial runs last and combines the < 100 tables into a single table. The size of this single table can be very large so you can optionally drop all columns but those you are directly interested in.

Again run as a batch job.
Reduce to the < 100 tables: slurm_merge_tables.sh_

.. _slurm_merge_tables.sh: https://github.com/legacysurvey/obiwan/blob/master/bin/slurm_merge_tables.sh

Merge the < 100 files into a single table: slurm_merge_tables_serial.sh_

.. _slurm_merge_tables_serial.sh: https://github.com/legacysurvey/obiwan/blob/master/bin/slurm_merge_tables_serial.sh>


Analysis and Plotting
----------------------

The majority of plots from Chp 4-5 of my thesis were made from the derived tables using this script: ``obiwan/py/obiwan/qa/plots_randomsprops_fluxdiff.py``

I'd recommend running on your laptop using one of the < 100 merged derived tables, since they are manageable size and are a random sampling of bricks so the plots you get should be representative. Once everything is working, run on the large single merged derived table from a NERSC login node.
