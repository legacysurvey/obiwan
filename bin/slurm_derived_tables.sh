#!/bin/bash -l

#SBATCH -p debug
#SBATCH -N 10
#SBATCH -t 00:10:00
#SBATCH --account=desi
#SBATCH -J randoms
#SBATCH -L SCRATCH,project
###SBATCH -C haswell

export doWhat=heatmap_table
#export doWhat=randoms_table
export dataset=dr3
#export dataset=dr5
export outdir=eboss_elg
#export outdir=elg_${dataset}_1000per
export thedate="02-19-2018"
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
    python $CSCRATCH/obiwan_code/obiwan/py/obiwan/runmanager/derived_tables.py \
    --doWhat ${doWhat} --dataset ${dataset} --data_dir $CSCRATCH/obiwan_out/${outdir} \
    --date ${thedate} \
    --bricks_fn ${bricks_fn} --nproc ${tasks} 

