"""Reads 1000s per brick tables and merges them

mpi4py cannot 'gather' lists of fits_tables, so each rank must write out
its own fits_table, then a later serial job must be used to merge the 
tables from each rank
"""

import numpy as np
import os
import sys
from glob import glob

try: 
    from astrometry.util.fits import fits_table, merge_tables
except ImportError:
    pass

def dir_for_mpi(derived_dir):
    return os.path.join(derived_dir,'merged_tmp')

def dir_for_serial(derived_dir):
    return os.path.join(derived_dir,'merged')

class RandomsTable(object):
    def __init__(self,derived_dir):
        self.derived_dir= derived_dir

class HeatmapTable(object):
    def __init__(self,derived_dir):
        self.derived_dir= derived_dir

    def run(self,bricks_to_merge):
        return self.datarelease(bricks_to_merge)
    
    def datarelease(self,bricks_to_merge):
        Tlist=[]
        for brick in bricks_to_merge:
            fn= os.path.join(self.derived_dir,brick[:3],brick,
                             'heatmap_datarelease.fits')
            #print('Reading %s' % fn)
            Tlist.append( fits_table(fn))
        return merge_tables(Tlist,columns='fillzero')

def main_mpi(bricks=[],doWhat=None,
             nproc=1,derived_dir=None):
    """
    Args:
        nproc: > 1 for mpi4py
        bricks: list of bricks
    """
    outfn= os.path.join(derived_dir,
                        'merged_%s.fits' % doWhat)
    if os.path.exists(outfn):
        print('Merged table already exists %s' % outfn)
        return

    if nproc > 1:
        from mpi4py.MPI import COMM_WORLD as comm
        bricks= np.array_split(bricks, comm.size)[comm.rank]
    else:
        class MyComm(object):
            def __init__(self):
                self.rank=0
                self.size=1
        comm= MyComm()

    if doWhat == 'randoms_table':
        tabMerger= RandomsTable(derived_dir)
    elif doWhat == 'heatmap_table':
        tabMerger= HeatmapTable(derived_dir)
    
    tab= tabMerger.run(bricks)
    tmpDir= dir_for_mpi(derived_dir) 
    try:
        os.makedirs(tmpDir)
    except OSError:
        pass
    outfn= os.path.join(tmpDir,'%s_rank%d.fits' % (doWhat,comm.rank))
    tab.writeto(outfn)
    print('Wrote %s' % outfn)

def main_serial(doWhat=None,derived_dir=None):
    """merges the rank tables that are stored in merge_tmp/"""
    saveDir= dir_for_serial(derived_dir) 
    try:
        os.makedirs(saveDir)
    except OSError:
        pass
    outfn= os.path.join(saveDir,'%s.fits' % doWhat)
    if os.path.exists(outfn):
        print('Merged table already exists %s' % outfn)
        return

    search=os.path.join(dir_for_mpi(derived_dir),
                        '%s_rank*.fits' % doWhat)
    tab_fns= glob(search)
    if len(tab_fns) == 0:
        raise ValueError('found nothing with search: %s' % search)
    tabs=[]
    for fn in tab_fns:
        tabs.append( fits_table(fn))
    print('Merging %d tables' % len(tabs))
    tab= merge_tables(tabs,columns='fillzero')
    tab.writeto(outfn)
    print('Wrote %s' % outfn)
    print('has %d rows' % len(tab))



if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--doWhat', type=str, choices=['randoms_table','heatmap_table'],required=True)
    parser.add_argument('--derived_dir', type=str, required=True) 
    parser.add_argument('--nproc', type=int, default=1, help='set to > 1 to run mpi4py') 
    parser.add_argument('--bricks_fn', type=str, default=None,
                        help='specify a fn listing bricks to run, or a single default brick will be ran') 
    parser.add_argument('--merge_rank_tables', action="store_true", default=False,help="set to merge the rank tables in the merge_tmp/ dir")
    args = parser.parse_args()
   
    if args.merge_rank_tables:
        main_serial(doWhat=args.doWhat,
                    derived_dir=args.derived_dir)
        sys.exit(0)

    # Bricks to run
    if args.bricks_fn is None:
        bricks= ['1266p292']
    else:
        bricks= np.loadtxt(args.bricks_fn,dtype=str)

    kwargs= vars(args)
    for dropCol in ['bricks_fn','merge_rank_tables']:
        del kwargs[dropCol]
    kwargs.update(bricks=bricks)

    main_mpi(**kwargs)
