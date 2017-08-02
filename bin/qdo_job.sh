#! /bin/bash

export brick="$1"
source ~/.bashrc_desiconda_orig

# Force MKL single-threaded
# https://software.intel.com/en-us/articles/using-threaded-intel-mkl-in-multi-thread-application
export MKL_NUM_THREADS=1

#fn=testobiwan_$brick.txt
#echo brick=$brick >> $fn


# Try limiting memory to avoid killing the whole MPI job...
#ulimit -a


#outdir=$SCRATCH/dr3more

#bri=$(echo $brick | head -c 3)
mkdir logs
log="logs/$brick.log"

#echo Logging to: $log
#echo Running on ${NERSC_HOST} $(hostname)

#echo -e "\nStarting on ${NERSC_HOST} $(hostname)\n" >> $log

python obiwan/kenobi.py -b ${brick} \
    -n 2 --DR 5 -o elg --add_sim_noise --zoom 1550 1650 1550 1650 \
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
