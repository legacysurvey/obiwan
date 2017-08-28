#!/bin/bash -l

#SBATCH -p debug
#SBATCH -N 1
#SBATCH -t 00:05:00
#SBATCH --account=desi
#SBATCH -J obiwan
#SBATCH -L SCRATCH,project
#SBATCH -C haswell

export brick=1238p245
export rowstart=0
#export do_skipids=no
#export db_table=obiwan_elg_9deg
export do_skipids=yes
export db_table=obiwan_9deg_doskip
export object=elg
export dataset=dr5
export nobj=300

usecores=4
threads=$usecores
#threads=1

# Load production env
source $CSCRATCH/obiwan_code/obiwan/bin/run_atnersc/prodenv_obiwan
export LEGACY_SURVEY_DIR=$obiwan_data/legacysurveydir_${dataset}

# Redirect logs
export bri=$(echo $brick | head -c 3)
export outdir=${obiwan_out}/${db_table}
export log=${outdir}/${object}/${bri}/${brick}/rs${rowstart}/log.${brick}
mkdir -p $(dirname $log)
echo Logging to: $log

# NERSC / Cray / Cori / Cori KNL things
export KMP_AFFINITY=disabled
export MPICH_GNI_FORK_MODE=FULLCOPY
export MKL_NUM_THREADS=1
#export OMP_NUM_THREADS=1
export OMP_NUM_THREADS=$threads

cd $obiwan_code/obiwan/py
export dataset=`echo $dataset | tr '[a-z]' '[A-Z]'`
srun -n 1 -c $usecores python obiwan/kenobi.py --dataset ${dataset} -b ${brick} \
                       --nobj ${nobj} --rowstart ${rowstart} -o ${object} \
                       --db_table ${db_table} \
                       --outdir $outdir --add_sim_noise \
                       --threads $threads  \
                       --do_skipids $do_skipids \
                       >> $log 2>&1

