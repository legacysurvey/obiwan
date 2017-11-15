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
export ra1=173.5
export ra2=176.5
export dec1=23.0
export dec2=26.0
export nrandoms=240000
export kdedir=$CSCRATCH/obiwan_out/randoms
export outdir=$CSCRATCH/obiwan_out/randoms/${name_for_run}

# NERSC / Cray / Cori / Cori KNL things
export KMP_AFFINITY=disabled
export MPICH_GNI_FORK_MODE=FULLCOPY
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

srun -n 32 -c 1 \
    python $obiwan_code/obiwan/py/obiwan/draw_radec_color_z.py \
    --dowhat sample --obj elg --ra1 ${ra1} --ra2 ${ra2} --dec1 ${dec1} --dec2 ${dec2} \
    --ndraws ${nrandoms} --kdedir ${kdedir} --outdir ${outdir}

