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

def derived_field_dir(brick,data_dir):
    return os.path.join(data_dir,'derived',
                        brick[:3],brick)

def uniform_obiwan_randoms(brick,data_dir):
    """Returns a uniform randoms and obiwan randoms table per brick"""
    search= os.path.join(data_dir,'tractor',
                         brick[:3],brick,
                         'rs*','tractor-%s.fits' % brick)
    rsdirs= glob(search)
    rsdirs= [os.path.dirname(dr)
             for dr in rsdirs]
    if len(rsdirs) == 0:
        raise ValueError('no rsdirs found: %s' % search)
    uniform,obi= [],[]
    for dr in rsdirs:
        simcat= fits_table((os.path.join(dr,'simcat-elg-%s.fits' % brick)
                            .replace('/tractor/','/obiwan/')))
        idsadded= fits_table((os.path.join(dr,'sim_ids_added.fits')
                            .replace('/tractor/','/obiwan/')))
        tractor= fits_table(os.path.join(dr,'tractor-%s.fits' % brick))
        # Uniform randoms
        assert(len(idsadded) == len(set(idsadded.id)))
        simcat.cut( pd.Series(simcat.id).isin(idsadded.id) )
        uniform.append(simcat)
        # Obiwan randoms
        tractor.cut(tractor.brick_primary)
        cols= np.array(tractor.get_columns())
        del_cols= cols[(pd.Series(cols)
                          .str.startswith('apflux_'))]
        for col in del_cols:
            tractor.delete_column(col)
        # nearest match in (ra2,dec2) for each point in (ra1,dec1)
        I,J,d = match_radec(simcat.ra,simcat.dec,
                            tractor.ra,tractor.dec, 1./3600,
                            nearest=True)
        assert(np.all(d <= 1./3600))
        tractor.cut(J)
        for simkey in ['id','ra','dec']:
            tractor.set('simcat_%s' % simkey,
                        simcat.get(simkey)[I])
        obi.append(tractor)
    return (merge_tables(uniform, columns='fillzero'),
            merge_tables(obi, columns='fillzero'))

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
            print('ccd_fns=',ccd_fns)
            for fn in ccd_fns:
                tmp=fits_table(fn)
                print('number=%d' % len(expids))
                expids=  np.hstack([expids,tmp.expid])
            expids= np.array(list(set(expids)))
            np.save(save_fn,expids)
            print('Wrote %s' % save_fn)

def final_expids(outdir):
    save_fn= os.path.join(outdir,'derived',
                          'ccds_unique.npy')
    if os.path.exists(save_fn):
        print('Skipping, already exists %s' % save_fn)
    else:
        bri_fns= glob(os.path.join(outdir,'derived',
                                    'ccds','ccds_*.npy'))
        expids=np.array([])
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
    parser.add_argument('--nproc', type=int, default=1, help='set to > 1 to run mpi4py')  
    args = parser.parse_args()

    #mpi_expids_per_bri(nproc=args.nproc,
    #                  data_dir=args.data_dir,outdir=args.outdir)
    final_expids(args.outdir)

