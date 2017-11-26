#!/bin/bash -l

#SBATCH -p debug
#SBATCH -N 1
#SBATCH -t 00:30:00
#SBATCH --account=desi
#SBATCH -J obiwan
#SBATCH -L SCRATCH,project
#SBATCH -C haswell

source $CSCRATCH/obiwan_code/obiwan/bin/run_atnersc/bashrc_obiwan
export name_for_run=elg_100deg2
export ra1=150.0
export ra2=160.0
export dec1=0.
export dec2=10.0
export nrandoms=2400000
export outdir=$CSCRATCH/obiwan_out/randoms/${name_for_run}

# NERSC / Cray / Cori / Cori KNL things
export KMP_AFFINITY=disabled
export MPICH_GNI_FORK_MODE=FULLCOPY
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

srun -n 32 -c 1 \
    python $obiwan_code/obiwan/py/obiwan/draw_radec_color_z.py \
    --obj elg --ra1 ${ra1} --ra2 ${ra2} --dec1 ${dec1} --dec2 ${dec2} \
    --ndraws ${nrandoms} --outdir ${outdir} --nproc 32

