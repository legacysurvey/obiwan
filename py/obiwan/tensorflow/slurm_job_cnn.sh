#!/bin/bash
#SBATCH -N 1
#SBATCH -C knl,quad,cache
#SBATCH -p debug
#SBATCH -J tf
#SBATCH -t 00:05:00

module load tensorflow/intel-head
#module load tensorflow/1.4.0rc0
export OMP_NUM_THREADS=68
export KMP_AFFINITY="granularity=fine,verbose,compact,1,0"
export KMP_SETTINGS=1
export KMP_BLOCKTIME=1
export isKNL=yes

date
srun -n 68 -c 1 --cpu_bind=cores python cnn.py
date
