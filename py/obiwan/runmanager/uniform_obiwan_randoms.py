"""
Using obiwan outputs, make a table of 'official' randoms per-brick
- uniform randoms: random ra,dec + geometry cut
- obiwan randoms: uniform randoms + recovered by tractor
"""

import numpy as np
import os
from glob import glob
import pandas as pd
 
from obiwan.common import stack_tables

try: 
    from astrometry.util.fits import fits_table
    from astrometry.libkd.spherematch import match_radec
except ImportError:
    pass

def derived_field_dir(brick,data_dir):
    return os.path.join(data_dir,'derived',
                        brick[:3],brick)

def uniform_randoms(brick,data_dir):
    search= os.path.join(data_dir,'obiwan',
                         brick[:3],brick,
                         'rs*','simcat-elg-%s.fits' % brick
    simcat_fns= glob(search)
    if len(simcat_fns) == 0:
        raise ValueError('no files found: %s' % search)
    idsadded_fns= [fn.replace('simcat-elg-%s' % brick,
                              'sim_ids_added') 
                   fn fn in simcat_fns]
    raise ValueError
    simcat= stack_tables(simcat_fns,textfile=False)
    idsadded= stack_tables(idsadded_fns,textfile=False)
    assert(len(idsadded) == len(set(idsadded.id)))
    simcat.cut( pd.Series(simcat.id).isin(idsadded.id) )
    return simcat

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
        fn= os.path.join(derived_field_dir(brick,data_dir),
                         'uniform_randoms.fits')
        if os.path.exists(fn):
            print('Skipping, already exists %s' % fn)
        else:
            uniform= uniform_randoms(brick,data_dir)
            try:
                os.makedirs(os.path.dirname(fn))
            except IOError:
                pass
            uniform.writeto(fn)
            print('Wrote %s' % fn)


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
        bricks= ['1211p060'] #['1126p220']
    else:
        bricks= np.loadtxt(args.bricks_fn,dtype=str)

    mpi_main(nproc=args.nproc,bricks=bricks,
             data_dir=args.data_dir)

