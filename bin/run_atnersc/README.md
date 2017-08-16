# Running Obiwan @ NERSC

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

### Setup
For all runs do
```sh
$ export obiwan_data=$CSCRATCH/obiwan_data
$ mkdir $obiwan_data
$ cd $obiwan_data
$ wget http://portal.nersc.gov/project/desi/users/kburleigh/obiwan/legacysurveydirs.tar.gz
$ tar -xzv legacysurveydirs.tar.gz
```
For eBOSS DR3 do
```sh
$ export LEGACY_SURVEY_DIR=${obiwan_data}legacysurveydir_ebossdr3
```
or instead for DR5
```sh
$ export LEGACY_SURVEY_DIR=${obiwan_data}legacysurveydir_dr5
```

### Git clone 2 repos
```sh
$ export obiwan_code=$CSCRATCH/obiwan_code
$ mkdir $obiwan_code
$ cd $obiwan_code
$ git clone https://github.com/legacysurvey/obiwan.git
$ git clone https://github.com/legacysurvey/legacypipe.git
$ cd legacypipe
```
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
You're now ready to do a production run! Edit this slurm job script appropriately
https://github.com/legacysurvey/obiwan/blob/master/bin/shifter_job.sh

Then submit it with
```sh
$ cd $obiwan_repo/obiwan/py
$ sbatch ../bin/shifter_job.sh
```
### Qdo (optional)
...fill in




