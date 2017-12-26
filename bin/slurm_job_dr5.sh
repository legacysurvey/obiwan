#!/bin/bash -l

#SBATCH -p debug
#SBATCH -N 1
#SBATCH -t 00:10:00
#SBATCH --account=desi
#SBATCH -J obiwan
#SBATCH -L SCRATCH,project
#SBATCH -C haswell

export name_for_run=elg_dr5
export randoms_db=obiwan_elg_dr5
export dataset=dr5
#export brick=1090p295
export brick=1091p297
export rowstart=0
export do_skipids=no
export do_more=no
export minid=1
export object=elg
export nobj=1000

usecores=4
threads=$usecores
#threads=1

# Load production env
source $CSCRATCH/obiwan_code/obiwan/bin/run_atnersc/bashrc_obiwan

set +x
export LEGACY_SURVEY_DIR=$obiwan_data/legacysurveydir_${dataset}

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
# Protect against astropy configs
#export XDG_CONFIG_HOME=/dev/shm
#srun -n $SLURM_JOB_NUM_NODES mkdir -p $XDG_CONFIG_HOME/astropy

cd $obiwan_code/obiwan/py
export dataset=`echo $dataset | tr '[a-z]' '[A-Z]'`
srun -n 1 -c $usecores python obiwan/kenobi.py --dataset ${dataset} -b ${brick} \
                       --nobj ${nobj} --rowstart ${rowstart} -o ${object} \
                       --randoms_db ${randoms_db} \
                       --outdir $outdir --add_sim_noise \
                       --threads $threads  \
                       --do_skipids $do_skipids \
                       --do_more $do_more --minid $minid \
                       >> $log 2>&1

