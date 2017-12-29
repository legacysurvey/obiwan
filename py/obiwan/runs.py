"""
Grabs, configures, and otherwise sets up a production run on eBOSS DR3 data or DR5 data
"""

import os
from glob import glob
import numpy as np

from obiwan.common import stack_tables
from astrometry.util.fits import fits_table, merge_tables

DOWNLOAD_ROOT = "http://portal.nersc.gov/project/desi/users/kburleigh/obiwan/"
NERSC_ROOT = DOWNLOAD_ROOT.replace("http://portal.nersc.gov/project/",
                                   "/global/project/projectdirs/")\
                          .replace("/users/","/www/users/")

MAX_CCDS=62

def in_eboss(T):
    # TODO: ccd corners instead of center
    x= np.logical_and.reduce((T.ra > 317.,T.ra < 360., T.dec > -2, T.dec < 2))
    y= np.logical_and.reduce((T.ra > 0.,T.ra < 45., T.dec > -5.,T.dec < 5.))
    z= np.logical_and.reduce((T.ra > 126., T.ra < 169.,T.dec > 14.,T.dec < 29.))
    return np.logical_or.reduce((x,y,z)) 

def add_str_arrays(lis):
    assert(len(lis) > 1)
    a= lis[0]
    for b in lis[1:]:
        a= np.char.add(a,b)
    return a
                  
def rm_duplicates(T):
    """survey-ccds fits_table"""
    T.pid= add_str_arrays([T.expnum.astype(str),
                           np.char.strip(T.image_filename)])
    keep=np.ones(len(T),bool)
    pids= set(T.pid)
    for cnt,pid in enumerate(pids):
        if (cnt+1) % 100 == 0:
            print('%d/%d' % (cnt+1,len(pids)))
        ind= np.where(T.pid == pid)[0]
        # More than 62 ccds have same expnum, must be dup
        if ind.size > MAX_CCDS:
            tmp_keep= np.zeros(len(T[ind]),bool)
            for ccdid in set(T[ind].ccdname):
                tmp_ind=  np.where(T[ind].ccdname == ccdid)[0]
                # duplicated ccdname
                if tmp_ind.size > 1:
                    tmp_keep[ tmp_ind[1:] ] =True # drop these
            keep[ ind[tmp_keep] ]= False # drop these by cut() method
    print('%d/%d are duplicates' % (np.where(keep == False)[0].size,len(keep)))
    return keep

if __name__ == "__main__":
    eboss_fns= glob(os.path.join(NERSC_ROOT,'configs/dr3eBOSS/additional_ccds',
                               "survey-ccds-*.fits.gz"))
    dr3_fns= glob(os.path.join(NERSC_ROOT,'configs/dr3eBOSS/dr3',
                             "survey-ccds-*.fits.gz"))
    assert(len(eboss_fns) > 0)
    assert(len(dr3_fns) > 0)

    dr3= stack_tables(dr3_fns,textfile=False)
    eboss= stack_tables(eboss_fns,textfile=False)

    dr3.set('pid', add_str_arrays([dr3.expnum.astype(str),
                                 np.char.strip(dr3.image_filename)]))
    eboss.set('pid', add_str_arrays([eboss.expnum.astype(str),
                                   np.char.strip(eboss.image_filename)]))
    dr3.cut( in_eboss(dr3))
    eboss.cut( in_eboss(eboss))

    T=  merge_tables([dr3,eboss], columns='fillzero')
    keep= rm_duplicates(T)
    T.cut(keep)
    name='survey-ccds-ebossDR3.fits'
    T.writeto(name)
    print('Wrote %s' % name)

    a=set(dr3.pid).union(set(eboss.pid))
    fn=[lin.split("decam/")[1] for lin in a]
    name='eboss_image_list.txt'
    with open(name,'w') as foo:
        for f in fn:
            foo.write('%s\n' % f)
    print('Wrote %s' % name) 

 
