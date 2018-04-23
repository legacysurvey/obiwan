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
from collections import defaultdict


from obiwan.runmanager.derived_tables import TargetSelection
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
            fn= self.table_fn(brick)
            if os.path.exists(fn):
                tab= fits_table(self.table_fn(brick))
                Tlist.append(tab)
            else:
                print('Skipping brick %s b/c randoms table doesnt exist' % fn)
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
        d=defaultdict(list)
        for brick in bricks_to_merge:
            if os.path.exists( self.table_fn(brick)):
                d['brickname'].append(brick)
                self.add_obiwan_summary(d,self.table_fn(brick),
                                        prefix='tractor_')
            else:
                print('Skipping brick %s b/c randoms table doesnt exist')
        # Save
        T=fits_table()
        T.set('brickname', np.array(d['brickname']).astype(np.string_))
        for key in ['n_injected','n_recovered',
                    'n_elg_ngc','n_elg_sgc']:
            T.set(key, np.array(d[key]).astype(np.int32))
        T.set('brick_area', np.array(d['brick_area']).astype(np.float32))
        for b in 'grz':
            T.set('galdepth_'+b, np.array(d['galdepth_'+b]).astype(np.float32))
        self.write_table(T,self.savefn)

    def add_obiwan_summary(self,summary_dict,randoms_fn,prefix=''):
        try:
            T = fits_table(randoms_fn)
        except OSError:
            raise OSError('could not open %s' % randoms_fn)
        
        isRec= T.obiwan_mask == 1
        TS= TargetSelection(prefix='tractor_') 
        is_elg_ngc= TS.run(T,'eboss_ngc')
        is_elg_sgc= TS.run(T,'eboss_sgc')
        summary_dict['n_injected'].append( len(T))
        summary_dict['n_recovered'].append( len(T[isRec]) )
        summary_dict['n_elg_ngc'].append( len(T[isRec & is_elg_ngc]) )
        summary_dict['n_elg_sgc'].append( len(T[isRec & is_elg_sgc]) )
        # FIXME: depends on the brick
        summary_dict['brick_area'].append( 0.25**2 )
        for band in 'grz':
            keep= np.isfinite(T.get(prefix+'galdepth_'+band))
            depth= np.median(T.get(prefix+'galdepth_'+band)[keep])
            summary_dict['galdepth_'+band].append( depth)
 
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

def randoms_subset_count_rsdirs_per_brick(rand_subset_fn):
    """for Hui, count number of rsdirs per brick to get which of the randoms subsets bricks are done/"""
    a=fits_table(rand_subset_fn)
    df=pd.DataFrame(dict(brick=pd.Series(a.unique_id).str.split("_").str[1],
                         rsdir=pd.Series(a.unique_id).str.split("_").str[2]))
    num_rsdirs= df.groupby(['brick']).agg(lambda x: len(set(x)))
    num_rsdirs= num_rsdirs.reset_index().sort_values(by='rsdir',ascending=False)
    fn= rand_subset_fn.replace('.fits','_count_rsdirs_per_brick.csv')
    num_rsdirs.to_csv(fn,index=False)
    print('Wrote %s' % fn)


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
    parser.add_argument('--count_rsdirs_per_brick', action="store_true", default=False,help="read existing randoms_subset table")
    args = parser.parse_args()
   
    if args.merge_rank_tables:
        kwargs= vars(args)
        for key in ['merge_rank_tables','nproc','bricks_fn',
                    'count_rsdirs_per_brick']:
            del kwargs[key]
        main_serial(**kwargs)
        sys.exit(0)

    if args.count_rsdirs_per_brick:
        fn= os.path.join(args.derived_dir,'merged','randoms_subset.fits')
        #fn= '/Users/kaylan1/Downloads/obiwan_plots/randoms_subset_10k.fits'
        randoms_subset_count_rsdirs_per_brick(fn)
        sys.exit(0)
   
    # Bricks to run
    if args.bricks_fn is None:
        bricks= ['1266p292']
    else:
        bricks= np.loadtxt(args.bricks_fn,dtype=str)

    kwargs= vars(args)
    for dropCol in ['bricks_fn','merge_rank_tables',
                    'randoms_subset','count_rsdirs_per_brick']:
        del kwargs[dropCol]
    kwargs.update(bricks=bricks)
    
    main_mpi(**kwargs)
