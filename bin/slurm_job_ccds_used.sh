#!/bin/bash -l

#SBATCH -q debug
#SBATCH -N 1
#SBATCH -t 00:10:00
#SBATCH --account=desi
#SBATCH -J dr3
#SBATCH -L SCRATCH,project

export name=elg_dr5_1000per
export data_dir=/global/project/projectdirs/cosmo/data/legacysurvey/dr5
export outdir=/global/cscratch1/sd/kaylanb/obiwan_out/${name}
export which=per_bri

# Load production env
source $CSCRATCH/obiwan_code/obiwan/bin/run_atnersc/bashrc_desiconda_new

# NERSC / Cray / Cori / Cori KNL things
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
# Protect against astropy configs
#export XDG_CONFIG_HOME=/dev/shm
#srun -n $SLURM_JOB_NUM_NODES mkdir -p $XDG_CONFIG_HOME/astropy

export coresper=1
let tasks=24*${SLURM_JOB_NUM_NODES}/${coresper}
srun -n ${tasks} -c ${coresper} \
    python -u $CSCRATCH/obiwan_code/obiwan/py/obiwan/qa/unique_ccds_dr3.py \
    --data_dir ${data_dir} --outdir ${outdir} \
    --which ${which} --nproc ${tasks}

