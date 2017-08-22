#! /bin/bash

# Example
# qdo launch obiwan 3 --cores_per_worker 4 --batchqueue debug --walltime 00:05:00 --script $obiwan_code/obiwan/bin/qdo_job_test.sh --keep_env

export brick="$1"
export rs="$2"

# Cannnot source $CSCRATCH/obiwan_code/obiwan/bin/run_atnersc/prodenv_obiwan 
# inside a qdo script
#################
# We need 3 directories for Obiwan
export obiwan_data=$CSCRATCH/obiwan_data
export obiwan_code=$CSCRATCH/obiwan_code
export obiwan_out=$CSCRATCH/obiwan_out

# Python environment, Ted Kisner's desiconda for imaging  
module use $obiwan_code/obiwan/etc/modulefiles
module load obiwan_conda
# Env vars for imaging pipeline 
for name in unwise_coadds unwise_coadds_timeresolved; do
    module load $name
done  
export DUST_DIR=$obiwan_data/dust

# Obiwan and imaging pipeline repos
export PYTHONPATH=$obiwan_code/theValidator:${PYTHONPATH}
export PYTHONPATH=$obiwan_code/legacypipe/py:${PYTHONPATH}
export PYTHONPATH=$obiwan_code/obiwan/py:${PYTHONPATH}
##################

# Choose dataset
export dataset=dr5
export LEGACY_SURVEY_DIR=$obiwan_data/legacysurveydir_${dataset}

echo brick,rs= $brick,$rs

cd $obiwan_code/obiwan/py/obiwan
export log=$obiwan_outdir/tests/log.${brick}_rs${rs}
mkdir $(dirname $log)
echo logging to $log
pytest --cov=test/ --cov-config=test/.coveragerc --ignore=test/end_to_end test \
    >> $log 2>&1
pytest --cov=test/end_to_end/test_datasets.py test/end_to_end/test_datasets.py \
    >> $log 2>&1

