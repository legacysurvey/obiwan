#! /bin/bash

# Example
# qdo launch obiwan 3 --cores_per_worker 4 --batchqueue debug --walltime 00:05:00 --script $obiwan_code/obiwan/bin/qdo_job_test.sh --keep_env

export brick="$1"
export rowstart="$2"
export object=elg
export dataset=dr5
export nobj=100

# Load env, env vars
source $CSCRATCH/obiwan_code/obiwan/bin/run_atnersc/prodenv_obiwan 
export LEGACY_SURVEY_DIR=$obiwan_data/legacysurveydir_${dataset}
# assert we have some new env vars
: ${obiwan_code:?}

# Redirect logs
export bri=$(echo $brick | head -c 3)
export outdir=${obiwan_out}/${bri}/$brick
export log=${outdir}/${object}/${bri}/${brick}/rs${rowstart}/log.${brick}
mkdir -p $(dirname $log)
echo Logging to: $log

# NERSC / Cray / Cori / Cori KNL things
export KMP_AFFINITY=disabled
export MPICH_GNI_FORK_MODE=FULLCOPY
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

cd $obiwan_code/obiwan/py
export dataset=`echo $dataset | tr '[a-z]' '[A-Z]'`
python obiwan/kenobi.py --dataset ${dataset} -b ${brick} \
                        --nobj ${nobj} --rowstart ${rowstart} -o ${object} --outdir $outdir \
                        --add_sim_noise  \
                        >> $log 2>&1



