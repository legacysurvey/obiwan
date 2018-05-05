#!/bin/bash -l

#SBATCH -p debug
#SBATCH -N 1
#SBATCH -t 00:20:00
#SBATCH --account=desi
#SBATCH -J randoms
#SBATCH -L SCRATCH,project
###SBATCH -C haswell

export doWhat=randoms
#export doWhat=summary
#export outdir=eboss_elg
#export outdir=elg_dr5_1000per
export outdir=elg_dr5_500per
#export outdir=cosmos_subsets/subset60
#export outdir=cosmos_subsets/subset64
#export outdir=cosmos_subsets/subset69
#export thedate="02-27-2018"
export thedate="03-05-2018"
export derived_dir="${CSCRATCH}/obiwan_out/${outdir}/derived_${thedate}"
export bricks_fn=${derived_dir}/bricks.txt
#export bricks_fn=${CSCRATCH}/obiwan_out/bricks100.txt

# Load production env
source $CSCRATCH/obiwan_code/obiwan/bin/run_atnersc/bashrc_desiconda_new

# NERSC / Cray / Cori / Cori KNL things
export KMP_AFFINITY=disabled
export MPICH_GNI_FORK_MODE=FULLCOPY
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

if [ "${NERSC_HOST}" = "edison" ]; then
    export num_cores=24
else
    export num_cores=32
fi

let tasks=${num_cores}*${SLURM_JOB_NUM_NODES}
srun -n 1 -c ${num_cores} \
    python -u $CSCRATCH/obiwan_code/obiwan/py/obiwan/runmanager/merge_tables.py \
    --doWhat ${doWhat} --derived_dir ${derived_dir} \
    --bricks_fn ${bricks_fn} --nproc ${tasks} \
    --merge_rank_tables --randoms_subset 

