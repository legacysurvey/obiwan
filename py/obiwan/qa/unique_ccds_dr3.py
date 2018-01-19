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

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--data_dir', type=str, required=True, 
                        help='path to obiwan/, tractor/ dirs') 
    parser.add_argument('--outdir', type=str, required=True, 
                        help='path to write output file') 
    parser.add_argument('--which', type=str, choices=['per_bri','final'],required=True) 
    parser.add_argument('--nproc', type=int, default=1, help='set to > 1 to run mpi4py')  
    args = parser.parse_args()

    if args.which == 'per_bri':
        mpi_expids_per_bri(nproc=args.nproc,
                          data_dir=args.data_dir,outdir=args.outdir)
    elif args.which == 'final':
        final_expids(args.outdir)

