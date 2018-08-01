********************
﻿Obiwan using Docker
********************

The obiwan docker images are stored on my dockerhub_. ``obiwan_rtd_jupter`` will run obiwan, a jupyter notebook for developing obiwan, and the make html command for building obiwan docs. If that one doesn't work for some reason then this one does everything but juypter ``obiwan_rtd``.

.. _dockerhub: https://hub.docker.com/r/kaylanb/desi/tags/


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

.. _jupyter-with-docker:

Jupyter notebook
----------------------------

You can run Jupiter notebook server using the docker image and use your laptop's browser to display the port. This way you can run any code you want without having any of it installed on your laptop but still use your laptop's browser for the notebook and the git repos on your laptop for development.::


  docker run -it -v /Users/kaylanb/repos_for_docker:/srv/repos_for_docker -p 8888:8888  kaylanb/desi:obiwan_rtd_jupyter /bin/bash
  cd /srv/repos_for_docker/obiwan/doc/nb/
  jupyter notebook --allow-root

Then copy and paste the link into your browser. Edit code from your browser and commit to GitHub from your laptops terminal, while not needing a single line of obiwan or legacypipe or tractor installed on your laptop!

.. _docs-with-docker:

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

.. _run-a-brick-with-docker:

Run a brick at NERSC
---------------------

The following shows how to run a single brick injecting 1000 sources into DECam imaging. Obiwan works for MzLS and BASS imaging as well, but to run a MzLS/BASS brick you'll need to set up a new `legacypipe_dir` as I done for DECaLS imaging. I built the Docker images uses these Dockerfiles_. At NERSC, the ``shifterimg`` command is the stand-in for ``docker``, which has some of the functionality of docker and allows docker images to be pulled onto the login nodes and mounted across a many node compute job.::

  ssh <user>@cori.nersc.gov
  shifterimg -v pull kaylanb/desi:obiwan_rtd_jupyter

.. _Dockerfiles: https://github.com/legacysurvey/obiwan/tree/master/dockerfiles

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
