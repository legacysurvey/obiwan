#!/bin/bash -l

#SBATCH -p debug
#SBATCH -N 10
#SBATCH -t 00:10:00
#SBATCH --account=desi
#SBATCH -J train
#SBATCH -L SCRATCH,project
#SBATCH -C haswell

export which=sim
let tasks=32*${SLURM_JOB_NUM_NODES}
export bricks_fn=${CSCRATCH}/obiwan_out/elg_dr5_coadds/partially_done_bricks.txt

# Load production env
source $CSCRATCH/obiwan_code/obiwan/bin/run_atnersc/bashrc_desiconda

# NERSC / Cray / Cori / Cori KNL things
export KMP_AFFINITY=disabled
export MPICH_GNI_FORK_MODE=FULLCOPY
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

srun -n ${tasks} -c 1 \
    python $CSCRATCH/obiwan_code/obiwan/py/obiwan/tensorflow/create_training.py \
    --which ${which} --bricks_fn ${bricks_fn} --nproc ${tasks}
