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

class TargetSelection(object):
    def __init__(self,eboss_or_desi,dataset):
        assert(eboss_or_desi in ['eboss','desi'])
        assert(dataset in DATASETS)
        self.eboss_or_desi= eboss_or_desi
        self.dataset= dataset

    def keep(self,tractor):
        if self.eboss_or_desi == 'eboss':
            return self.ebossIsElg(tractor)
        elif self.eboss_or_desi == 'desi':
            return self.desiIsElg(tractor)

    def desiIsElg(self,tractor):
        kw= dict(primary=tractor.brick_primary)
        for band,iband in [('g',1),('r',2),('z',4)]:
            if self.dataset == 'dr5':
                kw[band+'flux']= tractor.get('flux_'+band) / tractor.get('mw_transmission_'+band)
            elif self.dataset == 'dr3':
                kw[band+'flux']= tractor.decam_flux[:,iband] / tractor.decam_mw_transmission[:,iband]
        return self._desiIsElg(**kw)
    
    def _desiIsElg(self,gflux=None, rflux=None, zflux=None, 
                   primary=None):
        """VERBATIM from 
        https://github.com/desihub/desitarget/blob/master/py/desitarget/cuts.py
        
        Args:
            gflux, rflux, zflux, w1flux, w2flux: array_like
                The flux in nano-maggies of g, r, z, w1, and w2 bands.
            primary: array_like or None
                If given, the BRICK_PRIMARY column of the catalogue.
        Returns:
            mask : array_like. True if and only the object is an ELG
                target.
        """
        #----- Emission Line Galaxies
        if primary is None:
            primary = np.ones_like(gflux, dtype='?')
        elg = primary.copy()
        elg &= rflux > 10**((22.5-23.4)/2.5)                       # r<23.4
        elg &= zflux > rflux * 10**(0.3/2.5)                       # (r-z)>0.3
        elg &= zflux < rflux * 10**(1.6/2.5)                       # (r-z)<1.6

        # Clip to avoid warnings from negative numbers raised to fractional powers.
        rflux = rflux.clip(0)
        zflux = zflux.clip(0)
        elg &= rflux**2.15 < gflux * zflux**1.15 * 10**(-0.15/2.5) # (g-r)<1.15(r-z)-0.15
        elg &= zflux**1.2 < gflux * rflux**0.2 * 10**(1.6/2.5)     # (g-r)<1.6-1.2(r-z)

        return elg 
    
    def ebossIsElg(self,tractor):
        return keep

class RandomsTable(object):
    def __init__(self,derived_dir):
        self.derived_dir= derived_dir
    
    def run(self,bricks_to_merge,fitsfn):
        for eboss_or_desi in ['eboss']:
            for randoms_tab in ['uniform','obiwan_a','obiwan_b','obiwan_real']:
                merged= self.select_targets_and_merge(bricks_to_merge,
                                                      randoms_table=randoms_tab,
                                                      eboss_or_desi=eboss_or_desi)
                fn=fitsfn.replace('.fits','_%s_%s.fits' % (randoms_tab,eboss_or_desi))
                merged.writeto(fn)
                print('Wrote %s' % fn)

    def select_targets_and_merge(bricks_to_merge,randoms_table='',
                                 eboss_or_desi=''):
        TS= TargetSelection(eboss_or_desi)
        Tlist=[]
        for brick in bricks_to_merge:
            fn= os.path.join(self.derived_dir,brick[:3],brick,
                             'randoms_%s.fits' % randoms_table)
            tractor= fits_table(fn)
            tractor.cut(TS.keep(tractor))
            Tlist.append(tractor)
        return merge_tables(Tlist,columns='fillzero')




class HeatmapTable(object):
    def __init__(self,derived_dir):
        self.derived_dir= derived_dir

    def run(self,bricks_to_merge,savefn):
        return self.datarelease(bricks_to_merge,savefn)
    
    def datarelease(self,bricks_to_merge,fitsfn):
        Tlist=[]
        for brick in bricks_to_merge:
            fn= os.path.join(self.derived_dir,brick[:3],brick,
                             'heatmap_datarelease.fits')
            #print('Reading %s' % fn)
            Tlist.append( fits_table(fn))
        T= merge_tables(Tlist,columns='fillzero')
        T.writeto(fitsfn)
        print('Wrote %s' % fitsfn)

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

    if doWhat == 'randoms':
        tabMerger= RandomsTable(derived_dir)
    elif doWhat == 'heatmap':
        tabMerger= HeatmapTable(derived_dir)
    
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
