"""Reads 1000s per brick tables and merges them

mpi4py cannot 'gather' lists of fits_tables, so each rank must write out
its own fits_table, then a later serial job must be used to merge the 
tables from each rank
"""

import numpy as np
import os
import sys
from glob import glob
import pandas as pd

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
    def __init__(self,derived_dir,savefn,**kwargs):
        self.derived_dir= derived_dir
        self.savefn=savefn

    def run(self,bricks_to_merge):
        Tlist=[]
        for brick in bricks_to_merge:
            tab= fits_table(self.table_fn(brick))
            Tlist.append(tab)
        T= merge_tables(Tlist,columns='fillzero')
        T.writeto(self.savefn)
        print('Wrote %s' % self.savefn)

    def table_fn(self,brick):
        return os.path.join(self.derived_dir,brick[:3],brick,
                            'name.fits')

class RandomsTable(MergeTable):
    def __init__(self,derived_dir,savefn):
        super().__init__(derived_dir,savefn)

    def table_fn(self,brick):
        return os.path.join(self.derived_dir,brick[:3],brick,
                            'randoms.fits')


class SummaryTable(MergeTable):
    """In addition to merging over brick tables, compute avg quantities per brick
    
    derived table "randoms.fits" must exist. Joins the brick summary 
    quantities from a data release with a similar set from the 
    randoms.fits table. Each brick's table has one 
    row and all tables get merged to make the eatmap plots
    """
    def __init__(self,derived_dir,savefn):
        """
        Args:
            rank: mpi rank
        """
        super().__init__(derived_dir,savefn)

    def table_fn(self,brick):
        return os.path.join(self.derived_dir,brick[:3],brick,
                            'randoms.fits')

     def run(self,bricks_to_merge):
        Tlist=[]
        for brick in bricks_to_merge:
            tab= self.obiwan_summary(self.table_fn(brick),
                                     prefix='tractor_')
            Tlist.append(tab)
        T= merge_tables(Tlist,columns='fillzero')
        self.write_table(T,self.savefn)

    def obiwan_summary(self,randoms_fn,prefix=''):
        T = fits_table(randoms_fn)
        
        isRec= T.obiwan_mask == 1
        summary= defaultdict(list)
        summary['frac_recovered'].append( len(T[isRec])/ float(len(T)))
        for band in 'grz':
            summary['galdepth_'+band]= np.median(T.get(prefix+'galdepth_'+band))
        # Type
        T=fits_table()
        T.set('frac_recovered', np.array(summary['frac_recovered']).astype(np.float32))
        for b in 'grz':
            T.set('galdepth_'+b, np.array(summary[b+'galdepth']).astype(np.float32))
        return T 
 
    def write_table(self,tab,fn):
        if not os.path.exists(fn): 
            tab.writeto(fn)
            print('Wrote %s' % fn)


def main_mpi(doWhat=None,bricks=[],nproc=1,
             derived_dir=None):
    """
    Args:
        nproc: > 1 for mpi4py
        bricks: list of bricks
    """
    if nproc > 1:
        from mpi4py.MPI import COMM_WORLD as comm
        bricks= np.array_split(bricks, comm.size)[comm.rank]
    else:
        class MyComm(object):
            def __init__(self):
                self.rank=0
                self.size=1
        comm= MyComm()

    tmpDir= dir_for_mpi(derived_dir) 
    try:
        os.makedirs(tmpDir)
    except OSError:
        pass
 
    if doWhat == 'randoms':
        savefn= os.path.join(tmpDir,'randoms_rank%d.fits' % comm.rank)
        tabMerger= RandomsTable(derived_dir,savefn)
    elif doWhat == 'summary':
        savefn= os.path.join(tmpDir,'summary_rank%d.fits' % comm.rank)
        tabMerger= SummaryTable(derived_dir,savefn)
    tab= tabMerger.run(bricks)

def fits_table_cols(randoms_fn):
    columns= ['unique_id','ra','dec','obiwan_mask',
              'tractor_anymask_g','tractor_anymask_r','tractor_anymask_z',
              'tractor_flux_g','tractor_flux_r','tractor_flux_z',
              'tractor_mw_transmission_g','tractor_mw_transmission_r','tractor_mw_transmission_z',
              'tractor_psfdepth_g','tractor_psfdepth_r','tractor_psfdepth_z']
    T= fits_table(randoms_fn, columns=columns)
    T.set('brickname',(pd.Series(T.unique_id).str.split('_')
                                             .str[1].values
                                             .astype(str)))
    return T

def main_serial(doWhat=None,derived_dir=None,
                randoms_subset=False):
    """merges the rank tables that are stored in merge_tmp/"""
    saveDir= dir_for_serial(derived_dir) 
    try:
        os.makedirs(saveDir)
    except OSError:
        pass

    if doWhat == 'randoms':
        wild= "randoms_rank*.fits"
        if randoms_subset:
            outfn= os.path.join(saveDir,'randoms_subset.fits')
        else:
            outfn= os.path.join(saveDir,'randoms.fits')
    elif doWhat == 'summary':
        wild= "summary_rank*.fits"
        outfn= os.path.join(saveDir,"summary.fits")
    
    if os.path.exists(outfn):
        print('Merged table already exists %s' % outfn)
        return

    search=os.path.join(dir_for_mpi(derived_dir),
                        wild)
    tab_fns= glob(search)
    if len(tab_fns) == 0:
        raise ValueError('found nothing with search: %s' % search)
    tabs=[]
    for fn in tab_fns:
        if randoms_subset:
            T= fits_table_cols(fn)
        else:
            T= fits_table(fn)
        tabs.append(T)
    print('Merging %d tables' % len(tabs))
    tab= merge_tables(tabs,columns='fillzero')
    tab.writeto(outfn)
    print('Wrote %s' % outfn)
    print('has %d rows' % len(tab))

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--doWhat', type=str, choices=['randoms','summary'],required=True)
    parser.add_argument('--derived_dir', type=str, required=True) 
    parser.add_argument('--nproc', type=int, default=1, help='set to > 1 to run mpi4py') 
    parser.add_argument('--bricks_fn', type=str, default=None,
                        help='specify a fn listing bricks to run, or a single default brick will be ran') 
    parser.add_argument('--merge_rank_tables', action="store_true", default=False,help="set to merge the rank tables in the merge_tmp/ dir")
    parser.add_argument('--randoms_subset', action="store_true", default=False,help="make a merged table that is a subset of the randoms columns")
    args = parser.parse_args()
   
    if args.merge_rank_tables:
        kwargs= vars(args)
        for key in ['merge_rank_tables','nproc','bricks_fn']:
            del kwargs[key]
        main_serial(**kwargs)
        sys.exit(0)
   
    # Bricks to run
    if args.bricks_fn is None:
        bricks= ['1266p292']
    else:
        bricks= np.loadtxt(args.bricks_fn,dtype=str)

    kwargs= vars(args)
    for dropCol in ['bricks_fn','merge_rank_tables',
                    'randoms_subset']:
        del kwargs[dropCol]
    kwargs.update(bricks=bricks)
    
    main_mpi(**kwargs)
