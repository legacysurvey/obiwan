************************
Obiwan with Desiconda
************************

Create your conda environment
----------------------------------------------------------------

The conda environment I used to carry out the eBOSS ELG runs was made by doing the following::

  cd desiconda
  CONFIG=cori-gcc-py27 PREFIX=/global/cscratch1/sd/kaylanb/obiwan_desiconda_add_pytest make clean
  CONFIG=cori-gcc-py27 PREFIX=/global/cscratch1/sd/kaylanb/obiwan_desiconda_add_pytest make imaging
  # Use NERSC's No Machine software b/c this command takes a LONG Time
  ./install_imaging_cori-gcc-py27.sh 2>&1 | tee log_add_pytest


Git clone Obiwan and Download auxillary data
-------------------------------------------------

Git clone the obiwan and legacypipe repos to a directory "obiwan_code"::

  export obiwan_code=$CSCRATCH/obiwan_code
  mkdir $obiwan_code
  cd $obiwan_code
  git clone https://github.com/legacysurvey/obiwan.git

``obiwan`` works with the ``dr5_wobiwan`` branch of ``legacypipe``. This is because I needed to add < 5 lines to the DR5-version of legacypipe to get it to work with `obiwan`, but legacypipe master is constantly being updated so I just made my own branch in legacypipe that I new would be stable.::

  git clone https://github.com/legacysurvey/legacypipe.git
  cd legacypipe
  git fetch
  git checkout dr5_wobiwan

Wget dataset files::

  export obiwan_data=$CSCRATCH/obiwan_data
  mkdir $obiwan_data
  cd $obiwan_data
  wget http://portal.nersc.gov/project/desi/users/kburleigh/obiwan/legacysurveydirs.tar.gz
  tar -xzvf legacysurveydirs.tar.gz


Dust map files::

  mkdir -p $obiwan_data/dust/maps
  cd $obiwan_data/dust/maps
  wget -c http://portal.nersc.gov/project/cosmo/temp/dstn/travis-ci/maps/SFD_dust_4096_ngp.fits
  wget -c http://portal.nersc.gov/project/cosmo/temp/dstn/travis-ci/maps/SFD_dust_4096_sgp.fits

Now you are ready to set some environment vars. Some point to /project, some to ``desiproc``'s space, and others to the auxillary data you just downloaded::

  for name in legacysurvey unwise_coadds unwise_coadds_timeresolved dust;do
    module unload $name
  done

Note, my desiconda install has sense been purged on Cori CSCRATCH, but it *may* help to look at the following file which my bashrc.ext sourced and which setup my desiconda environment.
``https://github.com/legacysurvey/obiwan/blob/master/bin/run_atnersc/bashrc_obiwan``
