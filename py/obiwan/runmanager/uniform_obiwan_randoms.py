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
    search= os.path.join(data_dir,'obiwan',
                         brick[:3],brick,
                         'rs*')
    rsdirs= glob(search)
    if len(rsdirs) == 0:
        raise ValueError('no rsdirs found: %s' % search)
    uniform,obi= [],[]
    for dr in rsdirs:
        simcat= fits_table(os.path.join(dr,'simcat-elg-%s.fits' % brick))
        idsadded= fits_table(os.path.join(dr,'sim_ids_added.fits'))
        tractor= fits_table((os.path.join(dr,'tractor-%s.fits' % brick)
                             .replace('/obiwan/','/tractor/')))
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

def mpi_main(nproc=1,data_dir='./',
             bricks=[]):
    """

    Args:
        nproc: > 1 for mpi4py
        bricks: list of bricks
    """
    if nproc > 1:
        from mpi4py.MPI import COMM_WORLD as comm
        bricks= np.array_split(bricks, comm.size)[comm.rank]
    
    for brick in bricks:
        dr= derived_field_dir(brick,data_dir)
        try:
            os.makedirs(dr)
        except OSError:
            pass
        uniform_fn= os.path.join(dr,'uniform_randoms.fits')
        obi_fn= os.path.join(dr,'obiwan_randoms.fits')
        if ((os.path.exists(uniform_fn)) &
            (os.path.exists(obi_fn))):
            print('Skipping, already exists %s %s' % \
                  (uniform_fn,obi_fn))
        else:
            uniform,obi= uniform_obiwan_randoms(brick,data_dir)
            uniform.writeto(uniform_fn)
            print('Wrote %s' % uniform_fn)
            obi.writeto(obi_fn)
            print('Wrote %s' % obi_fn)


if __name__ == '__main__':
    #testcase_main()
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--data_dir', type=str, required=True, 
                        help='path to obiwan/, tractor/ dirs') 
    parser.add_argument('--nproc', type=int, default=1, help='set to > 1 to run mpi4py') 
    parser.add_argument('--bricks_fn', type=str, default=None, 
                        help='specify a fn listing bricks to run, or a single default brick will be ran') 
    args = parser.parse_args()
    
    # Bricks to run
    if args.bricks_fn is None:
        bricks= ['1266p292']
    else:
        bricks= np.loadtxt(args.bricks_fn,dtype=str)

    mpi_main(nproc=args.nproc,bricks=bricks,
             data_dir=args.data_dir)

