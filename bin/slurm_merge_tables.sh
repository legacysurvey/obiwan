#!/bin/bash -l

#SBATCH -p debug
#SBATCH -N 5
#SBATCH -t 00:10:00
#SBATCH --account=desi
#SBATCH -J randoms
#SBATCH -L SCRATCH,project
###SBATCH -C haswell

export doWhat=heatmap_table
#export doWhat=randoms_table
export outdir=eboss_elg
#export outdir=elg_dr5_1000per
export derived_dir="${CSCRATCH}/obiwan_out/${outdir}/derived_02-19-2018"
export bricks_fn=${CSCRATCH}/obiwan_out/${outdir}/bricks.txt
#export bricks_fn=${CSCRATCH}/obiwan_out/bricks100.txt

# Load production env
source $CSCRATCH/obiwan_code/obiwan/bin/run_atnersc/bashrc_desiconda_new

# NERSC / Cray / Cori / Cori KNL things
export KMP_AFFINITY=disabled
export MPICH_GNI_FORK_MODE=FULLCOPY
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

let tasks=32*${SLURM_JOB_NUM_NODES}
srun -n ${tasks} -c 1 \
    python $CSCRATCH/obiwan_code/obiwan/py/obiwan/runmanager/merge_tables.py \
    --doWhat ${doWhat} --derived_dir ${derived_dir} \
    --bricks_fn ${bricks_fn} --nproc ${tasks} 

