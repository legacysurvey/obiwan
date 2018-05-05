#!/bin/bash -l

#SBATCH -q debug
#SBATCH -N 10
#SBATCH -t 00:10:00
#SBATCH --account=desi
#SBATCH -J dr3
#SBATCH -L SCRATCH,project


export ccds_dir=/global/cscratch1/sd/kaylanb/obiwan_data
#export ccds_table=legacysurveydir_cosmos/survey-ccds-cosmos-subsets60-69.fits.gz
export ccds_table=legacysurveydir_dr3/survey-ccds-decals-used.fits.gz
#export ccds_table=legacysurveydir_dr3/survey-ccds-eboss-dr3plus-used.fits.gz
#export ccds_table=legacysurveydir_dr3/survey-ccds-extra-used.fits.gz
#export ccds_table=legacysurveydir_dr3/survey-ccds-nondecals-used.fits.gz

# Load production env
source $CSCRATCH/obiwan_code/obiwan/bin/run_atnersc/bashrc_desiconda_new_eboss

# NERSC / Cray / Cori / Cori KNL things
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
# Protect against astropy configs
#export XDG_CONFIG_HOME=/dev/shm
#srun -n $SLURM_JOB_NUM_NODES mkdir -p $XDG_CONFIG_HOME/astropy

let tasks=24*${SLURM_JOB_NUM_NODES}
srun -n ${tasks} -c 1 \
    python -u $CSCRATCH/obiwan_code/obiwan/py/obiwan/runmanager/update_imagefn_if_doesnt_exist.py \
    --ccds_table ${ccds_dir}/${ccds_table} --nproc ${tasks}

