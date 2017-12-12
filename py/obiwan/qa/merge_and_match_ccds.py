from glob import glob
import os
import numpy as np
import pandas as pd
import subprocess

from astrometry.libkd.spherematch import match_radec
from astrometry.util.fits import fits_table, merge_tables


def stilts(infits,outfits):
    STILTSCMD='java -jar /global/cscratch1/sd/kaylanb/eboss/stilts.jar '

    tmpstr = (STILTSCMD+' tmatch1 action=keep1 '+
                'matcher=exact values="pid"  '+
                'in='+infits+' ifmt=fits '+
                'out='+outfits+' ofmt=fits')
    print(tmpstr)
    subprocess.call(tmpstr, shell=True)

if __name__ == '__main__':
    outdir='/global/cscratch1/sd/kaylanb/eboss'
    
    fn='survey-ccds-dr3.fits'
    if not os.path.exists(os.path.join(outdir,fn)):
        dr3_dir= '/global/project/projectdirs/cosmo/data/legacysurvey/dr3'
        allccds= [fits_table(os.path.join(dr3_dir,
                                'survey-ccds-%s.fits.gz' % name))
                  for name in ['decals','extra','nondecals']]
        allccds= merge_tables(allccds, columns='fillzero')
        allccds.writeto(fn)
        print('Wrote %s' % fn)
    
    savefn='survey-ccds-ebossdr3-rank0.fits'
    if not os.path.exists(os.path.join(outdir,savefn)):
        fns= glob('/global/cscratch1/sd/raichoor/eBOSS_ELG_obiwan/legacysurvey-*-ccds.fits')
        assert(len(fns) == 13157)
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        fns= np.array_split(fns,comm.size)[comm.rank]
        allccds=[]
        for cnt,fn in enumerate(fns):
            if (cnt+1) % 10 == 0: 
                print('-----rank %d: %d/%d -----' % (comm.rank,cnt+1,len(fns)))
            allccds.append(fits_table(fn))
        allccds= merge_tables(allccds, columns='fillzero')
        savefn=os.path.join(outdir,savefn).replace('0.fits','%d.fits' % comm.rank)
        allccds.writeto(savefn)
        print('Rank %s: Wrote %s' % (comm.rank,savefn))

    savefn='survey-ccds-ebossdr3.fits'
    if not os.path.exists(os.path.join(outdir,savefn)):
        fns= glob(os.path.join(outdir,'survey-ccds-ebossdr3*.fits'))
        assert(len(fns) > 0)
        allccds=[]
        for cnt,fn in enumerate(fns):
            allccds.append(fits_table(fn))
        allccds= merge_tables(allccds, columns='fillzero')
        allccds.writeto(savefn)
        print('Wrote %s' % savefn)

    # pid
    savefn='survey-ccds-ebossdr3-pid.fits'
    if not os.path.exists(os.path.join(outdir,savefn)):
        eboss= fits_table(savefn.replace('-pid',''))
        expids= pd.Series(eboss.expid).str.split('-').str[0].values 
        ccdnames= pd.Series(eboss.ccdname).str.strip().values
        pid= (np.array(expids,dtype=object) +\
             np.array(ccdnames,dtype=object)).astype(str)
        a=fits_table()
        a.set('pid',pid)
        a.set('row',np.arange(len(eboss)))
        a.writeto(savefn)
        print('Wrote %s' % savefn)

    # Remove duplicates
    #savefn= os.path.join(outdir,str(os.getpid())+'.fits')
    #stilts(os.path.join(outdir,'survey-ccds-ebossdr3-pid.fits'),
    #      savefn)
    
    # match 
    pid= fits_table(os.path.join(outdir,'138559.fits'))
    eboss= fits_table(os.path.join(outdir,'survey-ccds-ebossdr3.fits'))
    print(len(eboss))
    eboss= eboss[pid.row]
    print(len(eboss))

    # 
