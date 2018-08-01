********************
﻿Obiwan using Docker
********************

Run obiwan with your laptop
----------------------------

Clone legacypipe and obiwan onto your laptop and edit, commit, push etc. from your laptop as usual.::

  mkdir ~/repos_for_docker; cd ~/repos_for_docker
  git clone https://github.com/legacysurvey/legacypipe.git
  cd legacypipe;git fetch;git checkout dr5_wobiwan
  cd ~/repos_for_docker
  git clone https://github.com/legacysurvey/obiwan.git


Now use the docker image to run the software while you edit, commit, push the code files from your laptop as usual. Install docker on your Mac, linux box, or windows, then::

  # WARNING: 7.5 GB image!
  docker pull kaylanb/desi:obiwan_rtd_jupyter
  docker run -v /Users/kaylanb/repos_for_docker:/srv/repos_for_docker -it kaylanb/desi:obiwan_rtd_jupyter /bin/bash
  # run the test suite
  python /srv/repos_for_docker/obiwan/tests/test_200x200_pixel_regions.py

Note, the directory ``/srv`` exists on the docker image but not ``/srv/repos_for_docker``. The latter is created when you use the -v option.

Versions of Tractor, Astrometry.net, and Legacypipe
"""""""""""""""""""""""""""""""""""""""""""""""""""""
I used ``tractor-dr5.2`` for my thesis work but the newest version of Tractor is installed on the docker image. Astrometry.net is installed on the docker image as  ``astrometry-0.72`` which is the version I used for my thesis.

Jupyter notebook
----------------------------

You can run Jupiter notebook server using the docker image and use your laptop's browser to display the port. This way you can run any code you want without having any of it installed on your laptop but still use your laptop's browser for the notebook and the git repos on your laptop for development.::


  docker run -it -v /Users/kaylanb/repos_for_docker:/srv/repos_for_docker -p 8888:8888  kaylanb/desi:obiwan_rtd_jupyter /bin/bash
  cd /srv/repos_for_docker/obiwan/doc/nb/
  jupyter notebook --allow-root

Then copy and paste the link into your browser. Edit code from your browser and commit to GitHub from your laptops terminal, while not needing a single line of obiwan or legacypipe or tractor installed on your laptop!

Docs
----------------------------

The Docker image handles all obiwan code execution including development  from your local machine and data exploration with Jupiter on your local machine's browser; however, ReadTheDocs does not support user created Docker images. So it is simplest for developers to build the docs themselves, commit the files to branch ``gh-pages``, push to GitHub, and then the docs will be displayed here:
`https://legacysurvey.github.io/obiwan`_

The one drawback is that there will not be separate docs for each branch. But oh well!

To build the docs and push them online, do::

  docker run -it -v /Users/kaylanb/repos_for_docker:/srv/repos_for_docker  kaylanb/desi:obiwan_rtd_jupyter /bin/bash
  cd /srv/repos_for_docker/obiwan/doc
  make html

The ``_build/html`` directory contains all the files the server needs for the docs webpage. Copy that entire directory to the ``gh-pages`` branch. On your laptop (NOT inside the docker image), do::

  mkdir ~/Downloads/docs
  cd ~/repos_for_docker/obiwan
  cp -r docs/_build/html/* ~/Downloads/docs/
  git checkout gh-pages
  cd ~/repos_for_docker/obiwan # must be in base repo directory
  cp -r ~/Downloads/docs/* ./
  git add -u :/
  git commit -m "updated docs"
  git push origin gh-pages

At NERSC
---------

I built the Docker images uses the Dockerfiles here. At NERSC, the ``shifterimg`` command is the stand-in for ``docker``, which has some of the functionality of docker and allows docker images to be pulled onto the login nodes and mounted across a many node compute job.::

  ssh <user>@cori.nersc.gov
  shifterimg -v pull kaylanb/desi:obiwan_rtd_jupyter

``shifter`` (no "img" suffix) is what runs the image as an executable, e.g. ``shifter`` is NERSC's version of ``docker run``.

Now let's setup a compute job. Git clone obiwan and legacyipe *into this specific directory*.::

  mkdir $CSCRATCH/repos_for_docker
  cd $CSCRATCH/repos_for_docker
  git clone https://github.com/legacysurvey/obiwan.git
  cd obiwan; git fetch; git checkout mzls_bass; cd ../
  git clone https://github.com/legacysurvey/legacypipe.git
  cd legacypipe;git fetch;git checkout dr5_wobiwan; cd ../

Then run your code from a difference directory using these batch jobs. Note, `job.sh` is necessary because the docker entry point file I used doesn't seem to be loaded by shifter, so do it manually instead of automatically when the docker image is called. The following will run the obiwan test suite using the docker image. the outputs are automatically written to directories in the obiwan repo named ``obiwan/tests/out_testcase*``::

  mkdir $CSCRATCH/docker_runs
  cd $CSCRATCH/docker_runs
  cp $CSCRATCH/repos_for_docker/obiwan/bin/docker_job_testcases* ./
  sbatch docker_job_testcases.slurm


Now lets do something real. Inject 1000 ELGs into brick 1351p192. The following will draw randoms and independently draw ELG properties from my eBOSS ELG properties distribution. The ELG property samples are in 4 csv files, one of which is 8MB so download them::

  cd $CSCRATCH/repos_for_docker/obiwan/etc
  wget http://portal.nersc.gov/project/desi/users/kburleigh/obiwan/eboss_elg_sample_csvs.tar.gz
  tar -xzf eboss_elg_sample_csvs.tar.gz
  rm eboss_elg_sample_csvs.tar.gz

Now draw the samples::

  # output dir
  mkdir $CSCRATCH/docker_runs/test_brick
  cd $CSCRATCH/docker_runs
  shifter --image=kaylanb/desi:obiwan_rtd_jupyter /bin/bash
  source /srv/py3_venv/bin/activate
  export PYTHONPATH=$CSCRATCH/repos_for_docker/obiwan/py:$CSCRATCH/repos_for_docker/legacypipe/py:$PYTHONPATH
  python $CSCRATCH/repos_for_docker/obiwan/py/obiwan/draw_radec_color_z.py --survey eboss --obj elg --ra1 135.0 --ra2 135.3 --dec1 19.1 --dec2 19.3 --ndraws 5000 --outdir test_brick

This writes a randoms table named "randoms_seed_1_startid_1.fits".

Now run a brick injecting ELGs from this randoms table::

  cd $CSCRATCH/docker_runs
  cp $CSCRATCH/repos_for_docker/obiwan/bin/docker_job_one_brick* ./
  sbatch docker_job_one_brick.slurm

Edit ``docker_job_one_brick.slurm`` to change the brick, output directory, randoms table etc. You should now have the following in the test_brick directory::

  test_brick
  |-- randoms
  |   |-- README.txt
  |   `-- randoms_seed_1_startid_1.fits
  |-- logs
  |   `-- 135
  |       `-- 1351p192
  |           `-- rs0
  |               `-- log.1351p192
  |-- obiwan
  |   `-- 135
  |       `-- 1351p192
  |           `-- rs0
  |               |-- metacat-elg-1351p192.fits
  |               |-- sim_ids_added.fits
  |               `-- simcat-elg-1351p192.fits
  |-- metrics
  |   `-- 135
  |       `-- 1351p192
  |           `-- rs0
  |               |-- all-models-1351p192.fits
  |               `-- blobs-1351p192.fits.gz
  |-- tractor
  |   `-- 135
  |       `-- 1351p192
  |           `-- rs0
  |               |-- brick-1351p192.sha256sum
  |               `-- tractor-1351p192.fits
  `-- tractor-i
  |    `-- 135
  |        `-- 1351p192
  |            `-- rs0
  |                `-- tractor-1351p192.fits
  |-- coadd
  |   `-- 135
  |       `-- 1351p192
  |           `-- rs0
  |               |-- legacysurvey-1351p192-ccds.fits
  |               |-- legacysurvey-1351p192-chi2-g.fits.fz
  |               |-- legacysurvey-1351p192-chi2-r.fits.fz
  |               |-- legacysurvey-1351p192-chi2-z.fits.fz
  |               |-- legacysurvey-1351p192-image-g.fits.fz
  |               |-- legacysurvey-1351p192-image-r.fits.fz
  |               |-- legacysurvey-1351p192-image-z.fits.fz
  |               |-- legacysurvey-1351p192-image.jpg
  |               |-- legacysurvey-1351p192-invvar-g.fits.fz
  |               |-- legacysurvey-1351p192-invvar-r.fits.fz
  |               |-- legacysurvey-1351p192-invvar-z.fits.fz
  |               |-- legacysurvey-1351p192-model-g.fits.fz
  |               |-- legacysurvey-1351p192-model-r.fits.fz
  |               |-- legacysurvey-1351p192-model-z.fits.fz
  |               |-- legacysurvey-1351p192-model.jpg
  |               |-- legacysurvey-1351p192-resid.jpg
  |               |-- legacysurvey-1351p192-sims-g.fits.fz
  |               |-- legacysurvey-1351p192-sims-r.fits.fz
  |               |-- legacysurvey-1351p192-sims-z.fits.fz
  |               `-- legacysurvey-1351p192-simscoadd.jpg

  25 directories, 31 files


Next Steps
-----------

The above shows how to run a single brick injecting N sources for DECam imaging. Obiwan works for MzLS and BASS imaging as well, but those legacypipe_dirs have not been setup yet. Also, Obiwan has not been updated to use the current HEAD of legacypipe; instead, it is using DR5-era legacypipe.

******************
Post-processing
******************

All post-processing of an obiwan production run is done by ``obiwan/py/obiwan/runmanager/derived_tables.py``
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

**********************
Running eBOSS ELGs
**********************

The point of this section is to summarize how I setup and completed the eBOSS ELG runs in Chp5 of my thesis. I'd checkout the thesis out for more info if needed.

The outputs from all of my runs are backed up on the tape archive (HPSS) at NERSC under my user, kaylanb. There must be a way to make these readable by other users.

- set up

 - link to existing how to pages

- run obiwan

  - use docker or the latest desiconda-imaging build

  - links to existing how to pages

- post-process

	 - ditto

- analysis

	 - ditto

Recommendations for finishing the eBOSS ELG run
------------------------------------------------

For finishing the bricks I didn't, I'd recommend first doing the following:

The module ``obiwan/db_tools.py`` has a function for querying the PSQL db for all randoms contained in a given brick, you should add a new function that takes the output of that query function and writes it to a FITS tabled named something like ``randoms/bri/brick.fits``. Then run obiwan using the ``--randoms_from_fits`` option to read the randoms to inject from a fits table instead of the PSQL db. This should be a pre-processing step. e.g. Make a list of all the bricks I didn't finish (that you want to finish), make this randoms table for each brick, then run obiwan on those bricks using those randoms tables. PSQL is an added complexity that we can (and should) avoid in the obiwan run.
