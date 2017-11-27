#!/bin/bash -l

#SBATCH -p debug
#SBATCH -N 1
#SBATCH -t 00:05:00
#SBATCH --account=desi
#SBATCH -J obiwan
#SBATCH -L SCRATCH,project
#SBATCH -C haswell

source $CSCRATCH/obiwan_code/obiwan/bin/run_atnersc/bashrc_obiwan
export name_for_run=elg_dr5
export ra1=109.0
export ra2=278.5
export dec1="-11.1"
export dec2=35.4
export nrandoms=32000000
#export nrandoms=0.032e9
#24 nodes
#export nrandoms=0.767e9
export outdir=$CSCRATCH/obiwan_out/randoms/${name_for_run}

# NERSC / Cray / Cori / Cori KNL things
export KMP_AFFINITY=disabled
export MPICH_GNI_FORK_MODE=FULLCOPY
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

let tasks=32*$SLURM_JOB_NUM_NODES
srun -n ${tasks} -c 1 \
    python $obiwan_code/obiwan/py/obiwan/draw_radec_color_z.py \
    --obj elg --ra1 ${ra1} --ra2 ${ra2} --dec1 "${dec1}" --dec2 ${dec2} \
    --ndraws ${nrandoms} --outdir ${outdir} --nproc ${tasks}

