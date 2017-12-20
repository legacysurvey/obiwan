# Background
[The Tractor](https://github.com/dstndstn/tractor) is the forward modeling code that finds the best point spread function (PSF) or galaxy model fit to astrophysical sources in imaging data. [Legacypipe](https://github.com/legacysurvey/legacypipe) is the our main pipeline which analyzes all the imaging data and that uses The Tractor. 

# Run at NERSC

Ted Kisner builds our software stack on Cori and Edison, such as The Tractor and Astrometry.net. You can use by adding the following to your bashrc
```sh
export NERSC_HOST=/usr/common/usg/bin/nersc_host
module use /global/common/cori/contrib/desi/desiconda/20170818-1.1.12-img/modulefiles
module load desiconda
```

The code knows about various data directories through environment variables. Set these
```sh
module use /global/cscratch1/sd/desiproc/modulefiles
module load legacypipe/legacysurvey
module load legacypipe/unwise_coadds
module load legacypipe/unwise_coadds_timeresolved
module load legacypipe/dust
```

Now git clone the Legacypipe repo
```sh
cd $CSCRATCH
git clone https://github.com/legacysurvey/legacypipe.git
cd legacypipe
git checkout tags/dr5.0
export PYTHONPATH=$CSCRATCH/legacypipe/py:$PYTHONPATH
```
and run the test cases
```
mkdir $CSCRATCH/legacypipe_out
cd $CSCRATCH/legacypipe_out
python $CSCRATCH/legacypipe/py/test/runbrick_test.py 
```

# Clone the existing desiconda install
you can also do 
```sh
module use /global/common/${NERSC_HOST}/contrib/desi/desiconda/20170818-1.1.12-spec/modulefiles
module load desiconda/20170818-1.1.12-spec
conda create --prefix $SCRATCH/desiconda_20170818-1.1.12-img --file $DESICONDA/pkg_list.txt
source activate /scratch2/scratchdirs/kaylanb/desiconda_20170818-1.1.12-img
```
