#!/bin/bash -l

#SBATCH -p debug
#SBATCH -N 1
#SBATCH -t 00:10:00
#SBATCH --account=desi
#SBATCH -J randoms
#SBATCH -L SCRATCH,project
###SBATCH -C haswell

#export doWhat=heatmap
export doWhat=targets
#export outdir=eboss_elg
#export outdir=elg_dr5_1000per
export outdir=elg_dr5_500per
export thedate="02-22-2018"
#export eboss_or_desi=eboss
export eboss_or_desi=desi
export dr_or_obiwan=datarelease
#export dr_or_obiwan=obiwan
export derived_dir="${CSCRATCH}/obiwan_out/${outdir}/derived_${thedate}"
export bricks_fn=${CSCRATCH}/obiwan_out/${outdir}/bricks.txt
#export bricks_fn=${CSCRATCH}/obiwan_out/bricks100.txt

# Load production env
source $CSCRATCH/obiwan_code/obiwan/bin/run_atnersc/bashrc_desiconda_new

# NERSC / Cray / Cori / Cori KNL things
export KMP_AFFINITY=disabled
export MPICH_GNI_FORK_MODE=FULLCOPY
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

for randoms_table in uniform obiwan_a obiwan_b obiwan_real; do 
    echo Running randoms_table ${randoms_table}
    srun -n 1 -c 1 \
        python $CSCRATCH/obiwan_code/obiwan/py/obiwan/runmanager/merge_tables.py \
        --doWhat ${doWhat} --derived_dir ${derived_dir} \
        --randoms_table ${randoms_table} --eboss_or_desi ${eboss_or_desi} \
        --dr_or_obiwan ${dr_or_obiwan} \
        --bricks_fn ${bricks_fn} --merge_rank_tables
done 

