# Running Obiwan @ NERSC

### Code and Data
Git clone the repos to a directory "obiwan_code"
```sh
export obiwan_code=$CSCRATCH/obiwan_code
mkdir $obiwan_code
cd $obiwan_code
git clone https://github.com/legacysurvey/obiwan.git
git clone https://github.com/legacysurvey/theValidator.git
git clone https://github.com/legacysurvey/legacypipe.git
cd legacypipe
git fetch
git checkout dr5_wobiwan 
```

Wget dataset files
```sh
export obiwan_data=$CSCRATCH/obiwan_data
mkdir $obiwan_data
cd $obiwan_data
wget http://portal.nersc.gov/project/desi/users/kburleigh/obiwan/legacysurveydirs.tar.gz
tar -xzvf legacysurveydirs.tar.gz
```

Dust map files
```sh
mkdir -p $obiwan_data/dust/maps
cd $obiwan_data/dust/maps
wget -c http://portal.nersc.gov/project/cosmo/temp/dstn/travis-ci/maps/SFD_dust_4096_ngp.fits
wget -c http://portal.nersc.gov/project/cosmo/temp/dstn/travis-ci/maps/SFD_dust_4096_sgp.fits
```

### Python environment
I created a conda environment for Obiwan using Ted Kisner's [desiconda](https://github.com/desihub/desiconda.git) package for the imaging pipeline. You activate it with a NERSC "module load" command, after which you have the usual conda functionality, like "conda info -e". I put the module load command and setting additional env vars in https://github.com/legacysurvey/obiwan/blob/master/bin/run_atnersc/bashrc_obiwan, so once you git clone obiwan you simply source this file. From a clean environment on Cori,
```
source $obiwan_code/obiwan/bin/run_atnersc/bashrc_obiwan
```
alternatively you can copy the file to your $HOME 
```sh
cp $obiwan_code/obiwan/bin/run_atnersc/bashrc_obiwan
```
then all you have to do is login to Cori and
```sh
source ~/bashrc_obiwan
```

### Run unit tests
It is a good idea to pull master before running anything
```sh
cd $obiwan_code/obiwan
git pull origin master
```

Now run the first suite of unit tests
```sh
cd $obiwan_code/obiwan/py/obiwan
pytest --cov=test/ --cov-config=test/.coveragerc --ignore=test/end_to_end test/
```
which should print "3 passed in ..." some number of seconds.


If the first suite works, run the second suite. 
```sh
pytest --cov=test/end_to_end/test_datasets.py test/end_to_end/test_datasets.py
```
which should print "2 passed" and possibly some "warnings" and "ERROR: Failed to generate report: No data to report." As long as you see "2 passed" it worked! 

The second suite does "end to end" tests, currently for datasets DR3 and DR5. Each injects 2 simulated ELGs into a 100x100 pixel cutout of a DECam CCD and runs the imaging pipeline on it. The outpurs are the usual imaging pipeline outputs (tractor catalogues, coadds, etc.) and obiwan metadata. 

If the second test suite printed "2 passed" then you can try some production runs.

### Login to the PostgreSQL DB
The psql db at NERSC is called "desi". It's hosted at scidb2.nersc.gov and you sign in with user "desi_user. You'll need the ".pgpass" config/password file (email Kaylan for it). Then put the ".pgass" file in $HOME/ on Cori. 

Make sure the bashrc_obiwan loaded the "postresql" module. Then do
```sh
psql -U desi_user -d desi -h scidb2.nersc.gov
```
It worked if that brings you to a "desi=>" prompt

### Production Runs

There are two basic ways of doing the production runs
 1) submit a slurm job script for every brick you want to run 
 2) automatically submit jobs with qdo

And two options telling the compute nodes about your python environment for obiwan
 * A) use the obiwan conda environment: source "prodenv_obiwan"
 * B) use a Docker image (coming soon)

#### 1A)
See https://github.com/legacysurvey/obiwan/blob/master/bin/slurm_job.sh

Edit these lines:
```sh
export brick=1238p245
export rowstart=0
export object=elg
export dataset=dr5
export nobj=100
```

Then run with
```sh
cd $obiwan_code
sbatch $obiwan_code/obiwan/bin/slurm_job.sh
```

#### 1B)
Coming soon

#### 2A)
See https://github.com/legacysurvey/obiwan/blob/master/bin/qdo_job.sh

Edit these lines:
```sh
export object=elg
export dataset=dr5
export nobj=100
```

Add list of bricks and indices of randoms as qdo tasks
```sh
cd $obiwan_code
for i in {100..500..100};do echo 1238p245 $i >> obiwan_qdo.txt;done
qdo load obiwan obiwan_qdo.txt
```

Now launch 5 qdo workers for the 5 qdo tasks you just made, using 6 hardware cores per task
```sh
cd $obiwan_code
qdo launch obiwan 5 --cores_per_worker 6 --batchqueue debug --walltime 00:30:00 --script $CSCRATCH/obiwan_code/obiwan/bin/qdo_job.sh --keep_env
```

#### 2B)
Coming soon

### Please ignore everything after this for now

The idea is for any NERSC user to easily do optimized production runs of Obiwan using a Docker Image. The steps are basically
 - unpack some tar.gz files to the user's scratch space
 - git clone https://github.com/kaylanb/obiwan.git
 - submit the included slurm job script that will load the Docker Image and run the obiwan repo


### Possible Production Runs
You can do the following runs with obiwan

| Run | Sources | Docs |
| ------ | ------ | ------ |
| eBOSS DR3 | ELG | fill in |
| DR5 | ELG | fill in |

### Notes
I made my conda environment by 
* cd desiconda
* CONFIG=cori-gcc-py27 PREFIX=/global/cscratch1/sd/kaylanb/obiwan_desiconda_add_pytest make clean
* CONFIG=cori-gcc-py27 PREFIX=/global/cscratch1/sd/kaylanb/obiwan_desiconda_add_pytest make imaging
* (from NX) ./install_imaging_cori-gcc-py27.sh 2>&1 | tee log_add_pytest

Legacypipe
For eBOSS dr3
```sh
$ git checkout tags/dr3e
```
or instead for dr5
For eBOSS dr3
```sh
$ git checkout tags/dr5.0
```

### If other desiconda modules already loaded:
```sh
for name in desiconda; do module unload $name;done
for name in legacysurvey unwise_coadds unwise_coadds_timeresolved dust;do 
  module unload $name
done
```

### Docker Image
Make sure you can see the correct Docker Image on NERSC Cori.You should see two images with this command
```sh
shifterimg images|grep tskisner|grep imaging
```
This one is for eBOSS DR3
```sh
cori       docker     READY    bcdab57daa   2017-08-04T19:50:54 tskisner/desiconda:1.1.9-imaging-py27
```
and this one is for DR5
```sh
cori       docker     READY    85235b9309   2017-08-03T13:49:09 tskisner/desiconda:1.1.9-imaging
```

### How I build the Docker images
Again using Ted Kisner's desiconda repo. First I installed docker-ce on ubuntu following instruction in `docker_install.sh`. Then using desiconda to make the Dockerfile. For py27,

```sh
cd desiconda_fork
git checkout add_pytest
CONFIG=docker-gcc-py27 PREFIX=$HOME/docker_images VERSION=20170822 make clean
CONFIG=docker-gcc-py27 PREFIX=$HOME/docker_images VERSION=20170822 make imaging
sudo ocker build --file Dockerfile_imaging_docker-gcc-py27 --tag py27:add-pytest . 2>&1 | tee log_py27_addpytest.txt 
```

### Run!
First we'll submit a single job and make sure it works. Later we'll install qdo to submit and manage 1000s of jobs. Simply do
```sh
cd $obiwan_code/obiwan/py
sbatch ../bin/run_atnersc/shifter_job.sh
```
There are two things you may want to edit in that script
 1) the #SBATCH lines
 2) the env var LEGACY_SURVEY_DIR, which you set to either "${obiwan_data}/legacysurveydir_ebossdr3" for eBOSS DR3 or to "${obiwan_data}/legacysurveydir_dr5" for DR5

Warning: make sure you submit the job from a clean environment

### Submit and Manage 1000s of jobs with Qdo
We'll use the conda package manager for python. We only want a basic environment to install qdo with, so Miniconda will suffice
```sh
export obiwan_conda=$CSCRATCH/obiwan_conda
mkdir $obiwan_conda
wget  https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash 
```
When it asks where to install, don't accept the default, instead choose
```sh
>>> /global/cscratch1/sd/<user>/obiwan_conda
```
and "Do you wish the installer to prepend the Miniconda3 install location
to PATH in your" bashrc, say no
```sh
>>> no
```
Activate your new environment and switch to python2.7 (obiwan is not python3 yet).
```sh
source $obiwan_conda/bin/activate
conda create -n py27 python=2.7 psycopg2
source activate py27
```
See the NERSC [docs](http://www.nersc.gov/users/data-analytics/data-analytics-2/python/anaconda-python/#toc-anchor-3). for more info. 

Now we can install qdo
```sh
cd $obiwan_repo
git clone https://bitbucket.org/berkeleylab/qdo.git
cd qdo
python setup.py install
cd ../
```
now type
```sh
qdo list
```
and you should see
```sh
QueueName              State   Waiting  Pending  Running Succeeded   Failed
edr                    Active       0        0        0       575        0
obiwan                 Active       0        0        0         5        1
dr5                    Active       0   123079     1760     51545      673
```

### Remaining
* setup desi_user db account 
* change obiwan to use desi_user not desi_admin for getSrcsInBrick

