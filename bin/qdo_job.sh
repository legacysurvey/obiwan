#! /bin/bash

export brick="$1"
source ~/.bashrc_desiconda_orig

# Force MKL single-threaded
# https://software.intel.com/en-us/articles/using-threaded-intel-mkl-in-multi-thread-application
export MKL_NUM_THREADS=1

echo brick=$brick
fn=env.txt
echo writing env vars to $fn
set > $fn
echo done

# Try limiting memory to avoid killing the whole MPI job...
#ulimit -a


#outdir=$SCRATCH/dr3more

#bri=$(echo $brick | head -c 3)
#mkdir -p $outdir/logs/$bri
#log="$outdir/logs/$bri/$brick.log"

#echo Logging to: $log
#echo Running on ${NERSC_HOST} $(hostname)

#echo -e "\nStarting on ${NERSC_HOST} $(hostname)\n" >> $log

#python legacypipe/runbrick.py \
#     --skip \
#     --threads 24 \
#     --skip-calibs \
#     --checkpoint checkpoints/${bri}/checkpoint-${brick}.pickle \
#     --pickle "pickles/${bri}/runbrick-%(brick)s-%%(stage)s.pickle" \
#     --brick $brick --outdir $outdir --nsigma 6 \
#     >> $log 2>&1

# qdo launch dr2n 16 --cores_per_worker 8 --walltime=24:00:00 --script ../bin/pipebrick.sh --batchqueue regular --verbose
# qdo launch edr0 4 --cores_per_worker 8 --batchqueue regular --walltime 4:00:00 --script ../bin/pipebrick.sh --keep_env --batchopts "--qos=premium -a 0-3"
