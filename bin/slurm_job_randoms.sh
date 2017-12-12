#!/bin/bash -l

#SBATCH -p debug
#SBATCH -N 1
#SBATCH -t 00:10:00
#SBATCH --account=desi
#SBATCH -J obiwan
#SBATCH -L SCRATCH,project
#SBATCH -C haswell

# USE PY3, py2 doesn't work with draw_points_eboss
source $CSCRATCH/obiwan_code/obiwan/bin/run_atnersc/bashrc_desiconda
export outdir=$CSCRATCH/obiwan_out/eboss_elg/randoms_test
export survey=eboss

export startid=1618001
export max_prev_seed=32
# SGC
# SGC A
#export nrandoms=418000
#export ra1=316.5
#export ra2=360.
#export dec1=-2.
#export dec2=2.
# SGC B
export nrandoms=1080000
export ra1=0.
export ra2=45.
export dec1=-5.
export dec2=5.
# NGC
#export nrandoms=1200000
#export ra1=126.
#export ra2=165.
#export dec1=14.
#export dec2=29.


# NERSC / Cray / Cori / Cori KNL things
export KMP_AFFINITY=disabled
export MPICH_GNI_FORK_MODE=FULLCOPY
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

#let tasks=32*$SLURM_JOB_NUM_NODES
tasks=16
srun -n ${tasks} -c 2 \
    python $obiwan_code/obiwan/py/obiwan/draw_radec_color_z.py \
    --survey ${survey} --obj elg --ra1 ${ra1} --ra2 ${ra2} --dec1 "-5" --dec2 ${dec2} \
    --ndraws ${nrandoms} --outdir ${outdir} --nproc ${tasks} \
    --startid ${startid} --max_prev_seed ${max_prev_seed}

