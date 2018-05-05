#!/bin/bash -l

#SBATCH -p debug
#SBATCH -N 1
#SBATCH -t 00:10:00
#SBATCH --account=desi
#SBATCH -J randoms
#SBATCH -L SCRATCH,project
#SBATCH -C haswell

# Load production env
source $CSCRATCH/obiwan_code/obiwan/bin/run_atnersc/bashrc_desiconda_new

# NERSC / Cray / Cori / Cori KNL things
export KMP_AFFINITY=disabled
export MPICH_GNI_FORK_MODE=FULLCOPY
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

let tasks=32*${SLURM_JOB_NUM_NODES}

srun -n ${tasks} -c 1 \
    python $CSCRATCH/obiwan_code/obiwan/py/obiwan/qa/unique_ccds_dr3.py \
    --data_dir /global/project/projectdirs/cosmo/data/legacysurvey/dr3 \
    --outdir /global/cscratch1/sd/kaylanb/obiwan_out/dr3_ccds_used \
    --which per_bri --nproc ${tasks}
