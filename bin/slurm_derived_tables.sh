#!/bin/bash -l

#SBATCH -p debug
#SBATCH -N 10
#SBATCH -t 00:10:00
#SBATCH --account=desi
#SBATCH -J randoms
#SBATCH -L SCRATCH,project
###SBATCH -C haswell

#export doWhat=heatmap
export doWhat=randoms
#export doWhat=targets
export eboss_or_desi=eboss
#export eboss_or_desi=desi
export db_randoms_table=obiwan_eboss_elg
#export db_randoms_table=obiwan_elg_dr5
export outdir=eboss_elg
#export outdir=elg_dr5_1000per
#export outdir=elg_dr5_500per
export thedate="02-27-2018"
export bricks_fn=${CSCRATCH}/obiwan_out/${outdir}/bricks.txt

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

srun -n ${tasks} -c 1 \
    python -u $CSCRATCH/obiwan_code/obiwan/py/obiwan/runmanager/derived_tables.py \
    --doWhat ${doWhat} --eboss_or_desi ${eboss_or_desi} --data_dir $CSCRATCH/obiwan_out/${outdir} \
    --db_randoms_table ${db_randoms_table} \
    --date ${thedate} \
    --bricks_fn ${bricks_fn} --nproc ${tasks} 

