#!/bin/bash -l

#SBATCH -p debug
#SBATCH -N 10
#SBATCH -t 00:10:00
#SBATCH --account=desi
#SBATCH -J randoms
#SBATCH -L SCRATCH,project
###SBATCH -C haswell

export doWhat=randoms
#export doWhat=heatmap
#export outdir=eboss_elg
#export outdir=elg_dr5_1000per
export outdir=elg_dr5_500per
#export outdir=cosmos_subsets/subset60
#export outdir=cosmos_subsets/subset64
#export outdir=cosmos_subsets/subset69
#export db_randoms_table=obiwan_eboss_elg
export db_randoms_table=obiwan_elg_dr5
#export db_randoms_table=obiwan_cosmos
#export dr3_or_dr5=dr3
export dr3_or_dr5=dr5
#export thedate="03-05-2018"
export thedate="03-16-2018"
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
    --doWhat ${doWhat} --dr3_or_dr5 ${dr3_or_dr5} --data_dir $CSCRATCH/obiwan_out/${outdir} \
    --db_randoms_table ${db_randoms_table} \
    --date ${thedate} \
    --bricks_fn ${bricks_fn} --nproc ${tasks} 

