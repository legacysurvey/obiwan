#!/bin/bash -l

#SBATCH -p debug
#SBATCH -N 1
#SBATCH -t 00:10:00
#SBATCH --account=desi
#SBATCH -J randoms
#SBATCH -L SCRATCH,project
###SBATCH -C haswell

#export doWhat=heatmap
export doWhat=randoms
#export doWhat=targets
export outdir=eboss_elg
#export outdir=elg_dr5_1000per
#export outdir=elg_dr5_500per
export thedate="02-27-2018"
export eboss_or_desi=eboss
#export eboss_or_desi=desi
#export dr_or_obiwan=datarelease
export dr_or_obiwan=obiwan
export derived_dir="${CSCRATCH}/obiwan_out/${outdir}/derived_${thedate}"

# Load production env
source $CSCRATCH/obiwan_code/obiwan/bin/run_atnersc/bashrc_desiconda_new

# NERSC / Cray / Cori / Cori KNL things
export KMP_AFFINITY=disabled
export MPICH_GNI_FORK_MODE=FULLCOPY
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

srun -n 1 -c 1 \
    python -u $CSCRATCH/obiwan_code/obiwan/py/obiwan/runmanager/merge_tables.py \
    --doWhat ${doWhat} --derived_dir ${derived_dir} \
    --eboss_or_desi ${eboss_or_desi} --dr_or_obiwan ${dr_or_obiwan} \
    --merge_rank_tables

