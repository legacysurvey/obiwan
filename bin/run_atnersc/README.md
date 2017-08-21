# Running Obiwan @ NERSC

### Setup
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
$ export obiwan_data=$CSCRATCH/obiwan_data
$ mkdir $obiwan_data
$ cd $obiwan_data
$ wget http://portal.nersc.gov/project/desi/users/kburleigh/obiwan/legacysurveydirs.tar.gz
$ tar -xzvf legacysurveydirs.tar.gz
```

Dust map files
```sh
$ mkdir -p $obiwan_data/dust/maps
$ cd $obiwan_data/dust/maps
$ wget -c http://portal.nersc.gov/project/cosmo/temp/dstn/travis-ci/maps/SFD_dust_4096_ngp.fits
$ wget -c http://portal.nersc.gov/project/cosmo/temp/dstn/travis-ci/maps/SFD_dust_4096_sgp.fits
```


Module load the conda environment for obiwan. I created it using Ted Kisner's [desiconda](https://github.com/desihub/desiconda.git) package and installing a few extra packages like "pytest". You simply grab the (bashrc_obiwan)[https://github.com/legacysurvey/obiwan/blob/master/bin/run_atnersc/bashrc_obiwan]
 file from the obiwan repo and source it. From a clean environment on Cori,
```
cd ~/
wget https://raw.githubusercontent.com/legacysurvey/obiwan/master/bin/run_atnersc/bashrc_obiwan
```
then whenever you want to run obiwan, ssh into a clean environment on Cori and
```sh
source ~/bashrc_obiwan
```

### Notes about the Obiwan bashrc file

You module load the obiwan conda environment with, 
```sh
module use $obiwan_code/obiwan/etc/modulefiles
module load obiwan_conda
```

These set the environment variables to run the leagcypipe pipeline
```sh
for name in dust unwise_coadds  unwise_coadds_timeresolved; do
    module load $name
done  
export LEGACY_SURVEY_DIR=$obiwan_data/legacysurveydir_dr5
```

this setups up your PYTHONPATH...
```sh
# Non-install pacakges
export PYTHONPATH=$obiwan_code/theValidator:${PYTHONPATH}
export PYTHONPATH=$obiwan_code/legacypipe/py:${PYTHONPATH}
export PYTHONPATH=$obiwan_code/obiwan/py:${PYTHONPATH}
```

### Run unit tests
This should return something like "3 passed in 5.44 seconds"
```sh
cd $obiwan_code/obiwan/py/obiwan
pytest --cov=test/ --cov-config=test/.coveragerc --ignore=test/end_to_end test/
```

Now run two full "end to end" tests, one using the DR3 dataset and the other for DR5. Each injects 2 simulated ELGs into a 100x100 pixel cutout of a DECam CCD, runs legacypipe on it, and writes out all the obiwan and legacypipe outputs. 
```sh
pytest --cov=test/end_to_end/test_datasets.py test/end_to_end/test_datasets.py
```
If the above worked you are ready for production runs!

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
* installing 201707...-imaging of desiconda, to 
`/global/cscratch1/sd/kaylanb/obiwan_desiconda/conda`
then
* pip install -U pytest
* pip install pytest-cov coveralls 

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

