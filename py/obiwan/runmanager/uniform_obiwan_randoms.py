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

DATASETS=['dr3','dr5']

def derived_field_dir(brick,data_dir,date):
    return os.path.join(data_dir,'derived_%s' % date,
                        brick[:3],brick)

def datarelease_dir(dataset):
    assert(dataset in DATASETS)
    proj='/global/project/projectdirs/cosmo/data/legacysurvey'
    return os.path.join(proj,dataset)
        

def obiwan_randoms_b(fn_obiwan_a,brick,dataset):
    """Computes one randoms table

    Args:
        fn_obiwan_a: obiwan_a randoms table fn for a given brick
        dataset: dr3,dr5 for the tractor cat of real sources
    
    Returns: 
        obiwan_b: obiwan_a but removing real sources, 
            e.g. sources with 1'' match in datarelease tractor cat
    """
    obiwan_a= fits_table(fn_obiwan_a)
    real= fits_table(os.path.join(datarelease_dir(dataset),
                                  'tractor',brick[:3],
                                  'tractor-%s.fits' % brick))
    # nearest match in (ra2,dec2) for each point in (ra1,dec1)
    I,J,d = match_radec(obiwan_a.ra,obiwan_a.dec,
                        real.ra,real.dec, 1./3600,
                        nearest=True)
    assert(np.all(d <= 1./3600))
    noMatch= np.ones(len(obiwan_a),dtype=bool)
    noMatch[I]= False
    obiwan_a.cut(noMatch)
    return obiwan_a

 

def uniform_obiwan_randoms(brick,data_dir):
    """Computes two randoms tables
    
    Returns: 
        uniform: random ra,dec cut to touching ccds
        obiwan_a: uniform that are recovered, e.g. having 1'' match in tractor cat 
    """
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

def mpi_main(nproc=1,data_dir='./',dataset=None,date='mm-dd-yyyy',
             bricks=[]):
    """

    Args:
        nproc: > 1 for mpi4py
        bricks: list of bricks
    """
    assert(dataset in DATASETS)
    if nproc > 1:
        from mpi4py.MPI import COMM_WORLD as comm
        bricks= np.array_split(bricks, comm.size)[comm.rank]
    
    for brick in bricks:
        dr= derived_field_dir(brick,data_dir,date)
        try:
            os.makedirs(dr)
        except OSError:
            pass
        fns= dict(uniform= os.path.join(dr,'uniform_randoms.fits'),
                  obiwan_a= os.path.join(dr,'obiwan_randoms_a.fits'),
                  obiwan_b= os.path.join(dr,'obiwan_randoms_b.fits'))
        if all((os.path.exists(fns[key]) 
                for key in fns.keys())):
            print('Skipping, already exist: ',fns)
        else:
            uniform,obiwan_a= uniform_obiwan_randoms(brick,data_dir)
            uniform.writeto(fns['uniform'])
            print('Wrote %s' % fns['uniform'])
            obiwan_a.writeto(fns['obiwan_a'])
            print('Wrote %s' % fns['obiwan_a'])
            obiwan_b= obiwan_randoms_b(fns['obiwan_a'],brick,dataset)
            obiwan_b.writeto(fns['obiwan_b'])
            print('Wrote %s' % fns['obiwan_b'])


if __name__ == '__main__':
    #testcase_main()
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--data_dir', type=str, required=True, 
                        help='path to obiwan/, tractor/ dirs') 
    parser.add_argument('--nproc', type=int, default=1, help='set to > 1 to run mpi4py') 
    parser.add_argument('--bricks_fn', type=str, default=None,
                        help='specify a fn listing bricks to run, or a single default brick will be ran') 
    parser.add_argument('--dataset', type=str, choices=['dr3','dr5'], 
                        help='for obiwan_randoms_b',required=True) 
    parser.add_argument('--date', type=str,help='mm-dd-yyyy, to label derived directory by',required=True) 
    args = parser.parse_args()
    
    # Bricks to run
    if args.bricks_fn is None:
        bricks= ['1266p292']
    else:
        bricks= np.loadtxt(args.bricks_fn,dtype=str)

    kwargs= dict(nproc=args.nproc,
                 bricks=bricks,
                 data_dir=args.data_dir,
                 dataset=args.dataset,
                 date=args.date)
    mpi_main(**kwargs)

