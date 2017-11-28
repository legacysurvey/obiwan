#! /bin/bash

# Example
# qdo launch obiwan 3 --cores_per_worker 4 --batchqueue debug --walltime 00:05:00 --script $obiwan_code/obiwan/bin/qdo_job_test.sh --keep_env

export name_for_run=elg_dr5_coadds
export randoms_db=obiwan_elg_dr5
export dataset=dr5
export brick="$1"
export rowstart="$2"
export do_skipids="$3"
export do_more="$4"
export object=elg
export nobj=300
export threads=4

# Load env, env vars
#source $CSCRATCH/obiwan_code/obiwan/bin/run_atnersc/prodenv_obiwan 
source $CSCRATCH/obiwan_code/obiwan/bin/run_atnersc/prodenv_desiconda_imaging
export LEGACY_SURVEY_DIR=$obiwan_data/legacysurveydir_${dataset}
# assert we have some new env vars
: ${obiwan_code:?}

# Redirect logs
export bri=$(echo $brick | head -c 3)
export outdir=${obiwan_out}/${name_for_run}
if [ ${do_skipids} == "no" ]; then
  if [ ${do_more} == "no" ]; then
    export rsdir=rs${rowstart}
  else
    export rsdir=more_rs${rowstart}
  fi
else
  if [ ${do_more} == "no" ]; then
    export rsdir=skip_rs${rowstart}
  else
    export rsdir=more_skip_rs${rowstart}
  fi
fi
export log=${outdir}/logs/${bri}/${brick}/${rsdir}/log.${brick}
mkdir -p $(dirname $log)
echo Logging to: $log

# NERSC / Cray / Cori / Cori KNL things
export KMP_AFFINITY=disabled
export MPICH_GNI_FORK_MODE=FULLCOPY
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

# Limit memory to avoid 1 srun killing whole node
if [ "$NERSC_HOST" = "edison" ]; then
    # 62 GB / Edison node = 65000000 kbytes
    maxmem=65000000
    let usemem=${maxmem}*${threads}/24
else
    # 128 GB / Cori node = 65000000 kbytes
    maxmem=134000000
    let usemem=${maxmem}*${threads}/32
fi
echo BEFORE
ulimit -Sa
ulimit -Sv $usemem
echo AFTER
ulimit -Sa


cd $obiwan_code/obiwan/py
export dataset=`echo $dataset | tr '[a-z]' '[A-Z]'`
python obiwan/kenobi.py --dataset ${dataset} -b ${brick} \
                        --nobj ${nobj} --rowstart ${rowstart} -o ${object} \
                        --randoms_db ${randoms_db} \
                        --outdir $outdir --add_sim_noise  \
                        --threads $threads  \
                        --do_skipids $do_skipids --do_more ${do_more} \
                        --early_coadds \
                        >> $log 2>&1



