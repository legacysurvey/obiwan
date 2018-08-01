#!/bin/bash -l

#SBATCH -A nstaff
#SBATCH -N 1 
#SBATCH -t 00:10:00
#SBATCH -q debug 
#SBATCH --image=kaylanb/desi:obiwan_rtd_jupyter
### Cori Haswell
#SBATCH -C haswell   
### Cori KNL
###SBATCH -C knl,quad,cache
### Cori
#SBATCH -L SCRATCH
### Edison
###SBATCH -L scratch3

date
shifter $CSCRATCH/docker_runs/job.sh
date

#export name_for_run=eboss_elg
#export randoms_db=obiwan_eboss_elg
#export dataset=dr3
#export brick=0304p045
#export rowstart=0
#export do_skipids=no
#export do_more=no
#export minid=1
#export object=elg
#export nobj=1000
#
#usecores=3
#threads=$usecores
##threads=1
#
## Load production env
#source $CSCRATCH/obiwan_code/obiwan/bin/run_atnersc/bashrc_desiconda_new_eboss
#
#export LEGACY_SURVEY_DIR=$obiwan_data/legacysurveydir_${dataset}
#
## Redirect logs
#export bri=$(echo $brick | head -c 3)
#export outdir=${obiwan_out}/${name_for_run}
#if [ ${do_skipids} == "no" ]; then
#  if [ ${do_more} == "no" ]; then
#    export rsdir=rs${rowstart}
#  else
#    export rsdir=more_rs${rowstart}
#  fi
#else
#  if [ ${do_more} == "no" ]; then
#    export rsdir=skip_rs${rowstart}
#  else
#    export rsdir=more_skip_rs${rowstart}
#  fi
#fi
#export log=${outdir}/logs/${bri}/${brick}/${rsdir}/log.${brick}
#mkdir -p $(dirname $log)
#echo Logging to: $log
#
## NERSC / Cray / Cori / Cori KNL things
#export MKL_NUM_THREADS=1
#export OMP_NUM_THREADS=1
## Protect against astropy configs
##export XDG_CONFIG_HOME=/dev/shm
##srun -n $SLURM_JOB_NUM_NODES mkdir -p $XDG_CONFIG_HOME/astropy
#
#cd $obiwan_code/obiwan/py
#srun -n 1 -c $usecores python obiwan/kenobi.py --dataset ${dataset} -b ${brick} \
#                       --nobj ${nobj} --rowstart ${rowstart} -o ${object} \
#                       --randoms_db ${randoms_db} \
#                       --outdir $outdir --add_sim_noise \
#                       --threads $threads  \
#                       --do_skipids $do_skipids \
#                       --do_more $do_more \
#                       --checkpoint \
#                       --skip_ccd_cuts \
#                       >> $log 2>&1

