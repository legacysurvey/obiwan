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

DATASETS=['dr3','dr5']

def dir_for_mpi(derived_dir):
    return os.path.join(derived_dir,'merged_tmp')

def dir_for_serial(derived_dir):
    return os.path.join(derived_dir,'merged')

class MergeTable(object):
    def __init__(self,derived_dir,**kwargs):
        self.derived_dir= derived_dir

    def run(self,bricks_to_merge,savefn):
        Tlist=[]
        for brick in bricks_to_merge:
            tab= fits_table(self.table_fn())
            Tlist.append(tab)
        T= merge_tables(Tlist,columns='fillzero')
        T.writeto(fitsfn)
        print('Wrote %s' % fitsfn)

    def table_fn(self):
        return os.path.join(self.derived_dir,brick[:3],brick,
                            'name.fits')
        

class TargetsTable(MergeTable):
    def __init__(self,derived_dir,
                 randoms_table,eboss_or_desi):
        super().__init__(derived_dir)
        assert(randoms_table in RANDOMS_TABLES)
        assert(eboss_or_desi in ['eboss','desi'])
        self.randoms_table= randoms_table
        self.eboss_or_desi= eboss_or_desi

    def table_fn(self):
        return os.path.join(self.derived_dir,brick[:3],brick,
                            'randoms_%s_%s.fits' % \
                            (self.randoms_table,self.eboss_or_desi))

class HeatmapTable(MergeTable):
    def __init__(self,derived_dir,dr_or_obiwan):
        super().__init__(derived_dir)
        assert(dr_or_obiwan in ['datarelease','obiwan'])
        self.dr_or_obiwan= dr_or_obiwan

    def table_fn(self):
        return os.path.join(self.derived_dir,brick[:3],brick,
                            'heatmap_%s.fits' % self.dr_or_obiwan)

def main_mpi(doWhat=None,bricks=[],nproc=1,
             derived_dir=None,
             randoms_table=None,eboss_or_desi=None,
             dr_or_obiwan=None)
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

    if doWhat == 'randoms':
        tabMerger= TargetsTable(derived_dir,
                                randoms_table,eboss_or_desi)
    elif doWhat == 'heatmap':
        tabMerger= HeatmapTable(derived_dir,
                                dr_or_obiwan)
    
    tmpDir= dir_for_mpi(derived_dir) 
    try:
        os.makedirs(tmpDir)
    except OSError:
        pass
    outfn= os.path.join(tmpDir,'%s_rank%d.fits' % (doWhat,comm.rank))
    tab= tabMerger.run(bricks,outfn)


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
    parser.add_argument('--doWhat', type=str, choices=['randoms','heatmap'],required=True)
    parser.add_argument('--derived_dir', type=str, required=True) 
    parser.add_argument('--randoms_table', type=str, default=None, choices=['uniform','obiwan_a',
                                            'obiwan_b','obiwan_real'],required=False)
    parser.add_argument('--eboss_or_desi', type=str, default=None, choices=['eboss','desi'],
                        required=False)
    parser.add_argument('--dr_or_obiwan', type=str, default=None, choices=['datarelease','obiwan'],
                        required=False)
    #parser.add_argument('--dataset', type=str, choices=['dr3','dr5'],required=True)
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
