#!/bin/bash -l
#SBATCH -M escori
#SBATCH -q xfer
#SBATCH -t 06:00:00
#SBATCH -J backup
#SBATCH -L SCRATCH

#squeue -M escori -u kaylanb 
export thedate="03_15_2018"
#export outdir=eboss_elg
#export outdir=elg_dr5_1000per
export outdir=elg_dr5_500per
#export outdir=cosmos_subsets/subset60
#export outdir=cosmos_subsets/subset64
#export outdir=cosmos_subsets/subset69

export refdir=$CSCRATCH/obiwan_out/${outdir}

cd $refdir
pwd
for name in obiwan randoms tractor tractor-i metrics logs checkpoint;do 
    echo $name
    htar -cf ${outdir}_${name}_${thedate}.tar ${name}
done

#cd coadd
#for dr in `ls`;do 
#    echo ${dr};htar -cf ${outdir}_coadd_${dr}_${thedate}.tar ${dr}
#done

