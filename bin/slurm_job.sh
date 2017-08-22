#!/bin/bash -l

#SBATCH -p debug
#SBATCH -N 1
#SBATCH -t 00:05:00
#SBATCH --account=desi
#SBATCH -J obiwan
#SBATCH -L SCRATCH,project
#SBATCH -C haswell

# Load production env
source $CSCRATCH/obiwan_code/obiwan/bin/run_atnersc/prodenv_obiwan
set -x
# Choose dataset
export dataset=dr5
export LEGACY_SURVEY_DIR=$obiwan_data/legacysurveydir_${dataset}
# Choose source to inject
export object=elg
# Choose row start index
export rowstart=0

export brick=1238p245
export bri=$(echo $brick | head -c 3)
export outdir=${obiwan_out}/${bri}/$brick
mkdir -p $outdir

usecores=4
threads=$usecores
#threads=1

export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=$threads

year=`date|awk '{print $NF}'`
today=`date|awk '{print $3}'`
month=`date +"%F"|awk -F "-" '{print $2}'`
export log=$outdir/${object}/${bri}/${brick}/rs${rowstart}/log.${brick}
mkdir -p $(dirname $log)
echo Logging to: $log

export dataset=`echo $dataset | tr '[a-z]' '[A-Z]'`
cd $obiwan_code/obiwan/py
srun -n 1 -c $usecores python obiwan/kenobi.py --dataset ${dataset} -b ${brick} \
                       -n 2 --rowstart ${rowstart} -o ${object} --outdir $outdir \
                       --add_sim_noise  \
                       >> $log 2>&1

echo SLURM FINISHED
