#!/bin/bash -l

#SBATCH -p debug
#SBATCH -N 10
#SBATCH -t 00:10:00
#SBATCH --account=desi
#SBATCH -J randoms
#SBATCH -L SCRATCH,project
#SBATCH -C haswell

export outdir=elg_dr5_500per
let tasks=32*${SLURM_JOB_NUM_NODES}
export bricks_fn=${CSCRATCH}/obiwan_out/${outdir}/bricks.txt

# Load production env
source $CSCRATCH/obiwan_code/obiwan/bin/run_atnersc/bashrc_desiconda_new

# NERSC / Cray / Cori / Cori KNL things
export KMP_AFFINITY=disabled
export MPICH_GNI_FORK_MODE=FULLCOPY
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

srun -n ${tasks} -c 1 \
    python $CSCRATCH/obiwan_code/obiwan/py/obiwan/runmanager/uniform_obiwan_randoms.py \
    --data_dir $CSCRATCH/obiwan_out/${outdir} --bricks_fn ${bricks_fn} --nproc ${tasks}
