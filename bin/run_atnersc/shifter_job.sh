#!/bin/bash
#SBATCH --image=docker:tskisner/desiconda:1.1.9-imaging-py27
#SBATCH --nodes=1
#SBATCH --partition=debug
#SBATCH -t 00:05:00
#SBATCH --account=desi
#SBATCH -J docker
#SBATCH -L SCRATCH,project,projecta
#SBATCH -C haswell

export brick=1238p245

# paths to pipeline code
export obiwan_code=$CSCRATCH/obiwan_code
export PYTHONPATH=$obiwan_code/obiwan/py:${PYTHONPATH}
export PYTHONPATH=$obiwan_code/legacypipe/py:${PYTHONPATH}
export PYTHONPATH=$obiwan_code/theValidator:${PYTHONPATH}
export PYTHONPATH=.:${PYTHONPATH}

# env vars for legacypipe pipeline
module use $obiwan_code/obiwan/etc/modulefiles
for name in dust unwise_coadds  unwise_coadds_timeresolved; do
    module load $name
done  

# images to use
export obiwan_data=$CSCRATCH/obiwan_data
#export LEGACY_SURVEY_DIR=${obiwan_data}/legacysurveydir_ebossdr3
export LEGACY_SURVEY_DIR=${obiwan_data}/legacysurveydir_dr5

# output directory for results
export obiwan_outdir=$CSCRATCH/obiwan_outdir
log="${obiwan_outdir}/logs/${brick}.log"
mkdir -p $(dirname $log)

# Force MKL single-threaded
# https://software.intel.com/en-us/articles/using-threaded-intel-mkl-in-multi-thread-application
export MKL_NUM_THREADS=1

echo logging to $log
srun -n 1 -c 32 shifter python obiwan/kenobi.py -b ${brick} \
    -n 2 --DR 5 -o elg --outdir $obiwan_outdir \
    --add_sim_noise --zoom 1550 1650 1550 1650 \
    >> $log 2>&1

#python legacypipe/runbrick.py \
#     --skip \
#     --threads 24 \
#     --skip-calibs \
#     --checkpoint checkpoints/${bri}/checkpoint-${brick}.pickle \
#     --pickle "pickles/${bri}/runbrick-%(brick)s-%%(stage)s.pickle" \
#     --brick $brick --outdir $outdir --nsigma 6 \
#     >> $log 2>&1

# qdo launch obiwan 3 --cores_per_worker 10 --batchqueue debug --walltime 00:05:00 --script $CSCRATCH/test/obiwan/bin/qdo_job.sh --keep_env
# qdo launch edr0 4 --cores_per_worker 8 --batchqueue regular --walltime 4:00:00 --script ../bin/pipebrick.sh --keep_env --batchopts "--qos=premium -a 0-3"
