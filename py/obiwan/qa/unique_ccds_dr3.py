"""
Using obiwan outputs, make a table of 'official' randoms per-brick
- uniform randoms: random ra,dec + geometry cut
- obiwan randoms: uniform randoms + recovered by tractor
"""

import numpy as np
import os
from glob import glob
import pandas as pd

try: 
    from astrometry.util.fits import fits_table, merge_tables
    from astrometry.libkd.spherematch import match_radec
except ImportError:
    pass

def mpi_expids_per_bri(nproc=1,data_dir='./',outdir='./'):
    """

    Args:
        nproc: > 1 for mpi4py
        bricks: list of bricks
    """
    bri_dirs= glob(os.path.join(data_dir,'coadd','*'))
    if nproc > 1:
        from mpi4py.MPI import COMM_WORLD as comm
        bri_dirs= np.array_split(bri_dirs, comm.size)[comm.rank]
    
    for dr in bri_dirs:
        new_dr= os.path.join(outdir,'derived','ccds')
        try:
            os.makedirs(new_dr)
        except OSError:
            pass
        save_fn= os.path.join(new_dr,'ccds_%s.npy' % os.path.basename(dr))
        if os.path.exists(save_fn):
            print('Skipping, already exists %s' % save_fn)
        else:
            expids=np.array([])
            ccd_fns= glob(os.path.join(dr,'*/legacysurvey-*-ccds.fits'))
            #print('ccd_fns=',ccd_fns)
            for fn in ccd_fns:
                try:
                    tmp=fits_table(fn)
                    if not 'expid' in tmp.get_columns():
                        tmp.set('expid',np.array(['%d-%s' % (a.expnum,a.ccdname) for a in tmp]))
                    #print('number=%d' % len(expids))
                    expids=  np.hstack([expids,tmp.expid])
                except OSError:
                    print('brick ccd corrupted: %s' % fn)
                    with open('ccds_used_corrupted.txt','a') as foo:
                        foo.write('brick ccd corrupted: %s\n' % fn)
            expids= np.array(list(set(expids)))
            np.save(save_fn,expids)
            print('Wrote %s' % save_fn)

def final_expids(outdir):
    save_fn= os.path.join(outdir,'derived',
                          'ccds_unique.npy')
    if os.path.exists(save_fn):
        print('Skipping, already exists %s' % save_fn)
    else:
        search= os.path.join(outdir,'derived',
                             'ccds','ccds_*.npy')
        bri_fns= glob(search)
        expids=np.array([])
        if len(bri_fns) == 0:
            raise ValueError('search returned nothing %s' % search)
        print('found %d bri_fns' % len(bri_fns))
        for fn in bri_fns:
            expids= np.hstack([expids,np.load(fn)])
        expids= np.array(list(set(expids)))
        np.save(save_fn,expids)
        print('Wrote %s' % save_fn)

def surveyccds_cut_to_expids(expid_npy_fn,survey_ccds_fn):
    """
    Args:
        expid_npy_fn: derived/ccds_unique.npy
        survey_ccds_fn: legacysurveydir_dr5/survey-ccds-legacypipe_wgain.fits.gz
    """
    expid_used=np.load(expid_npy_fn)
    ccds= fits_table(survey_ccds_fn,columns=['ra','dec','expnum','ccdname'])
    print('looping for expid')
    expid= np.array(['%s-%s' % (ccd.expnum,ccd.ccdname) for ccd in ccds])
    keep= pd.Series(np.char.strip(expid)).isin(np.char.strip(expid_used))
    ccds.cut(keep)
    savefn= expid_npy_fn.replace('.npy','.fits')
    ccds.writeto(savefn)
    print('Wrote %s' % savefn)


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--data_dir', type=str, required=True, 
                        help='path to obiwan/, tractor/ dirs') 
    parser.add_argument('--outdir', type=str, required=True, 
                        help='path to write output file') 
    parser.add_argument('--which', type=str, choices=['per_bri','final','survey_ccds'],
                        required=True) 
    parser.add_argument('--nproc', type=int, default=1, help='set to > 1 to run mpi4py')  
    parser.add_argument('--expid_npy_fn', required=False,help='if which == survey_ccds') 
    parser.add_argument('--survey_ccds_fn', required=False,help='if which == survey_ccds') 
    args = parser.parse_args()

    if args.which == 'per_bri':
        mpi_expids_per_bri(nproc=args.nproc,
                           data_dir=args.data_dir,outdir=args.outdir)
    elif args.which == 'final':
        final_expids(args.outdir)

    elif args.which == 'survey_ccds':
        surveyccds_cut_to_expids(args.expid_npy_fn,args.survey_ccds_fn) 

