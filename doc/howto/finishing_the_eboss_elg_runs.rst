************************************************
Finishing the eBOSS ELG Runs (ch5 of my thesis)
************************************************

Data is on HPSS
----------------

The outputs from all of my eBOSS ELG runs are backed up on the NERSC tape archive (HPSS) for user ``kaylanb``. I wrote the outputs to Cori CSCRATCH but nearly all of them have been purged. You can look here to see what's left ``/global/cscratch1/sd/kaylanb/obiwan_out]ls eboss_elg/``. There must be a way to make these readable by other users. In addition the derived tables are there. All the plots can be remade using the derived tables.

- outputs on tape: named ``eboss_elg*may7.tar``
- derived table: named ``eboss_elg_derived_04_11_2018_merged.tar``

You can find out which bricks need to be finished by looking at the derived tables. Once you finish them you'll need to start the **wonderful** task of transferring the terrabyotes of my finished brick outputs off tape, untar them, and merge with your terrabytes of outputs. Then make a final derived table that includes yours + my finished brick outputs.

Joint Distributions for eBOSS ELG flux, shape, and redshift
------------------------------------------------------------

Samples of eBOSS ELG flux, shape, and redshift are stored in 4 ``.csv`` files. Obiwan creates the randoms tables using `draw_radec_color_z.py` and it needs these ``.csv`` files. Running the obiwan test suite should tell you where to download them from. But just in case, there are here:
``http://portal.nersc.gov/project/desi/users/kburleigh/obiwan/eboss_elg_sample_csvs.tar.gz``
They are also on my google drive in case data on /project gets deleted for some reason.

Development that should be done before running the unfinished bricks
--------------------------------------------------------------------------------

The module ``obiwan/db_tools.py`` has a function for querying the PSQL db for all randoms contained in a given brick, you should add a new function that takes the output of that query function and writes it to a FITS tabled named something like ``randoms/bri/brick.fits``. Then run obiwan using the ``--randoms_from_fits`` option to read the randoms to inject from a fits table instead of the PSQL db. This should be a pre-processing step. e.g. Make a list of all the bricks I didn't finish (that you want to finish), make this randoms table for each brick, then run obiwan on those bricks using those randoms tables. PSQL is an added complexity that we can (and should) avoid in the obiwan run.
