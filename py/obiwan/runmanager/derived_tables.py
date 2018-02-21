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
        
class RandomsTables(object):
    """Creates the uniform,obiwan_a,obiwan_b randoms tables for a single brick"""
    def __init__(self, data_dir,dataset,date='mm-dd-yyyy'):
        assert(dataset in DATASETS)
        self.data_dir= data_dir
        self.dataset= dataset
        self.date= date

    def uniform_obiwana_obiwanb(self,brick):
        derived_dir= derived_field_dir(brick,self.data_dir,self.date):
        fns= dict(uniform= os.path.join(derived_dir,'uniform_randoms.fits'),
                  obiwan_a= os.path.join(derived_dir,'obiwan_randoms_a.fits'),
                  obiwan_b= os.path.join(derived_dir,'obiwan_randoms_b.fits'))
        if all((os.path.exists(fns[key]) 
                for key in fns.keys())):
            print('Skipping, already exist: ',fns)
        else:
            uniform,obiwan_a= self.uniform_obiwan_randoms(brick,self.data_dir)
            uniform.writeto(fns['uniform'])
            print('Wrote %s' % fns['uniform'])
            obiwan_a.writeto(fns['obiwan_a'])
            print('Wrote %s' % fns['obiwan_a'])
            obiwan_b= self.obiwan_randoms_b(fns['obiwan_a'],brick,self.dataset)
            obiwan_b.writeto(fns['obiwan_b'])
            print('Wrote %s' % fns['obiwan_b'])

    def uniform_obiwana(self,brick,data_dir):
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

    def obiwanb(self,fn_obiwan_a,brick,dataset):
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


class HeatmapTable(object):
    """Create a fits table with the heatmap values, one per brick"""
    def __init__(self, data_dir,dataset,date='mm-dd-yyyy'):
        assert(dataset in DATASETS)
        self.data_dir= data_dir
        self.dataset= dataset
        self.date= date

    def write_table_for_datarelease(self,brick):
        derived_dir= derived_field_dir(brick,self.data_dir,self.date):
        fn= os.path.join(derived_dir,'heatmap_datarelease.fits')
        if os.path.exists(fn): 
            print('Skipping, already exist: ',fn)
        else:
            T= self.get_table_for_datarelease([brick],self.dataset)
            T.writeto(fn)
            print('Wrote %s' % fn)
    
    def get_table_for_datarelease(self,bricklist,dataset):
        """
        Args:
            bricklist: Give a single brick as a list of length 1, e.g. [brick]
        """
        brickset = set()
        gn = []
        rn = []
        zn = []
        
        gnhist = []
        rnhist = []
        znhist = []
        
        nnhist = 6
        
        gdepth = []
        rdepth = []
        zdepth = []
        
        ibricks = []
        nsrcs = []
        npsf  = []
        nsimp = []
        nrex = []
        nexp  = []
        ndev  = []
        ncomp = []

        gpsfsize = []
        rpsfsize = []
        zpsfsize = []

        gpsfdepth = []
        rpsfdepth = []
        zpsfdepth = []
        ggaldepth = []
        rgaldepth = []
        zgaldepth = []

        wise_nobs = []
        wise_trans = []
        
        ebv = []
        gtrans = []
        rtrans = []
        ztrans = []
        
        
        #sfd = SFDMap()
        
        W = H = 3600
        # H=3600
        # xx,yy = np.meshgrid(np.arange(W), np.arange(H))
        unique = np.ones((H,W), bool)
        tlast = 0
       
        dirprefix= datarelease_dir(dataset)
        for ibrick,brick in enumerate(bricklist):
            #words = fn.split('/')
            #dirprefix = '/'.join(words[:-4])
            #print('Directory prefix:', dirprefix)
            #words = words[-4:]
            #brick = words[2]
            #print('Brick', brick)
            try:
                tfn = os.path.join(dirprefix, 'tractor', brick[:3], 'tractor-%s.fits'%brick)
                print('Tractor filename', tfn)
                if dataset == 'dr5':
                    T = fits_table(tfn, columns=['brick_primary', 'type',
                                                 'psfsize_g', 'psfsize_r', 'psfsize_z',
                                                 'psfdepth_g', 'psfdepth_r', 'psfdepth_z',
                                                 'galdepth_g', 'galdepth_r', 'galdepth_z',
                                                 'ebv',
                                                 'mw_transmission_g', 'mw_transmission_r', 'mw_transmission_z',
                                                 'nobs_w1', 'nobs_w2', 'nobs_w3', 'nobs_w4',
                                                 'mw_transmission_w1', 'mw_transmission_w2', 'mw_transmission_w3', 'mw_transmission_w4'])
                else:
                    T = fits_table(tfn, columns=['brick_primary', 'type', 'decam_psfsize',
                                             'decam_depth', 'decam_galdepth',
                                             'ebv', 'decam_mw_transmission',
                                             'wise_nobs', 'wise_mw_transmission'])
            except:
                print('Failed to read FITS table', tfn)
                import traceback
                traceback.print_exc()
                print('Carrying on.')
                continue

            brickset.add(brick)
            bricklist.append(brick)
            gn.append(0)
            rn.append(0)
            zn.append(0)

            gnhist.append([0 for i in range(nnhist)])
            rnhist.append([0 for i in range(nnhist)])
            znhist.append([0 for i in range(nnhist)])

            index = -1
            ibrick = np.nonzero(bricks.brickname == brick)[0][0]
            ibricks.append(ibrick)

            T.cut(T.brick_primary)
            nsrcs.append(len(T))
            types = Counter([t.strip() for t in T.type])
            npsf.append(types['PSF'])
            nsimp.append(types['SIMP'])
            nrex.append(types['REX'])
            nexp.append(types['EXP'])
            ndev.append(types['DEV'])
            ncomp.append(types['COMP'])
            print('N sources', nsrcs[-1])

            if dataset == 'dr5':
                gpsfsize.append(np.median(T.psfsize_g))
                rpsfsize.append(np.median(T.psfsize_r))
                zpsfsize.append(np.median(T.psfsize_z))

                gpsfdepth.append(np.median(T.psfdepth_g))
                rpsfdepth.append(np.median(T.psfdepth_r))
                zpsfdepth.append(np.median(T.psfdepth_z))

                ggaldepth.append(np.median(T.galdepth_g))
                rgaldepth.append(np.median(T.galdepth_r))
                zgaldepth.append(np.median(T.galdepth_z))

                wise_nobs.append(np.median(
                    np.vstack((T.nobs_w1, T.nobs_w2, T.nobs_w3, T.nobs_w4)).T,
                    axis=0))
                wise_trans.append(np.median(
                    np.vstack((T.mw_transmission_w1,
                               T.mw_transmission_w2,
                               T.mw_transmission_w3,
                               T.mw_transmission_w4)).T,
                               axis=0))

                gtrans.append(np.median(T.mw_transmission_g))
                rtrans.append(np.median(T.mw_transmission_r))
                ztrans.append(np.median(T.mw_transmission_z))
                
            else:
                gpsfsize.append(np.median(T.decam_psfsize[:,1]))
                rpsfsize.append(np.median(T.decam_psfsize[:,2]))
                zpsfsize.append(np.median(T.decam_psfsize[:,4]))

                gpsfdepth.append(np.median(T.decam_depth[:,1]))
                rpsfdepth.append(np.median(T.decam_depth[:,2]))
                zpsfdepth.append(np.median(T.decam_depth[:,4]))

                ggaldepth.append(np.median(T.decam_galdepth[:,1]))
                rgaldepth.append(np.median(T.decam_galdepth[:,2]))
                zgaldepth.append(np.median(T.decam_galdepth[:,4]))

                wise_nobs.append(np.median(T.wise_nobs, axis=0))
                wise_trans.append(np.median(T.wise_mw_transmission, axis=0))

                gtrans.append(np.median(T.decam_mw_transmission[:,1]))
                rtrans.append(np.median(T.decam_mw_transmission[:,2]))
                ztrans.append(np.median(T.decam_mw_transmission[:,4]))
                
            ebv.append(np.median(T.ebv))

            br = bricks[ibrick]

            #print('Computing unique brick pixels...')
            pixscale = 0.262/3600.
            wcs = Tan(br.ra, br.dec, W/2.+0.5, H/2.+0.5,
                      -pixscale, 0., 0., pixscale,
                      float(W), float(H))
            unique[:,:] = True
            self.find_unique_pixels(wcs, W, H, unique,
                                    br.ra1, br.ra2, br.dec1, br.dec2)
            U = np.flatnonzero(unique)
            #print(len(U), 'of', W*H, 'pixels are unique to this brick')
             
            index = bricklist.index(brick)
            assert(index == len(bricklist)-1)
        
            filepart = words[-1]
            filepart = filepart.replace('.fits.gz', '')
            filepart = filepart.replace('.fits.fz', '')
            print('File:', filepart)
            band = filepart[-1]
            assert(band in 'grz')
        
            nlist,nhist = dict(g=(gn,gnhist), r=(rn,rnhist), z=(zn,znhist))[band]
        
            upix = fitsio.read(fn).flat[U]
            med = np.median(upix)
            print('Band', band, ': Median', med)
            nlist[index] = med
        
            hist = nhist[index]
            for i in range(nnhist):
                if i < nnhist-1:
                    hist[i] = np.sum(upix == i)
                else:
                    hist[i] = np.sum(upix >= i)
            assert(sum(hist) == len(upix))
            print('Number of exposures histogram:', hist)
        
        ibricks = np.array(ibricks)
        
        #print('Maximum number of sources:', max(nsrcs))
        
        T = fits_table()
        T.brickname = np.array(bricklist)
        T.ra  = bricks.ra [ibricks]
        T.dec = bricks.dec[ibricks]
        T.nexp_g = np.array(gn).astype(np.int16)
        T.nexp_r = np.array(rn).astype(np.int16)
        T.nexp_z = np.array(zn).astype(np.int16)
        T.nexphist_g = np.array(gnhist).astype(np.int32)
        T.nexphist_r = np.array(rnhist).astype(np.int32)
        T.nexphist_z = np.array(znhist).astype(np.int32)
        T.nobjs  = np.array(nsrcs).astype(np.int16)
        T.npsf   = np.array(npsf ).astype(np.int16)
        T.nsimp  = np.array(nsimp).astype(np.int16)
        T.nrex   = np.array(nrex ).astype(np.int16)
        T.nexp   = np.array(nexp ).astype(np.int16)
        T.ndev   = np.array(ndev ).astype(np.int16)
        T.ncomp  = np.array(ncomp).astype(np.int16)
        T.psfsize_g = np.array(gpsfsize).astype(np.float32)
        T.psfsize_r = np.array(rpsfsize).astype(np.float32)
        T.psfsize_z = np.array(zpsfsize).astype(np.float32)
        with np.errstate(divide='ignore'):
            T.psfdepth_g = (-2.5*(-9.+np.log10(5.*np.sqrt(1. / np.array(gpsfdepth))))).astype(np.float32)
            T.psfdepth_r = (-2.5*(-9.+np.log10(5.*np.sqrt(1. / np.array(rpsfdepth))))).astype(np.float32)
            T.psfdepth_z = (-2.5*(-9.+np.log10(5.*np.sqrt(1. / np.array(zpsfdepth))))).astype(np.float32)
            T.galdepth_g = (-2.5*(-9.+np.log10(5.*np.sqrt(1. / np.array(ggaldepth))))).astype(np.float32)
            T.galdepth_r = (-2.5*(-9.+np.log10(5.*np.sqrt(1. / np.array(rgaldepth))))).astype(np.float32)
            T.galdepth_z = (-2.5*(-9.+np.log10(5.*np.sqrt(1. / np.array(zgaldepth))))).astype(np.float32)
        for k in ['psfdepth_g', 'psfdepth_r', 'psfdepth_z', 'galdepth_g', 'galdepth_r', 'galdepth_z']:
            v = T.get(k)
            v[np.logical_not(np.isfinite(v))] = 0.
        T.ebv = np.array(ebv).astype(np.float32)
        T.trans_g = np.array(gtrans).astype(np.float32)
        T.trans_r = np.array(rtrans).astype(np.float32)
        T.trans_z = np.array(ztrans).astype(np.float32)
        T.ext_g = -2.5 * np.log10(T.trans_g)
        T.ext_r = -2.5 * np.log10(T.trans_r)
        T.ext_z = -2.5 * np.log10(T.trans_z)
        T.wise_nobs = np.array(wise_nobs).astype(np.int16)
        T.trans_wise = np.array(wise_trans).astype(np.float32)
        T.ext_w1 = -2.5 * np.log10(T.trans_wise[:,0])
        T.ext_w2 = -2.5 * np.log10(T.trans_wise[:,1])
        T.ext_w3 = -2.5 * np.log10(T.trans_wise[:,2])
        T.ext_w4 = -2.5 * np.log10(T.trans_wise[:,3])
        return T 

    def find_unique_pixels(self,wcs, W, H, unique, ra1,ra2,dec1,dec2):
        if unique is None:
            unique = np.ones((H,W), bool)
        # scan the outer annulus of pixels, and shrink in until all pixels
        # are unique.
        step = 10
        for i in range(0, W//2, step):
            nu,ntot = self._ring_unique(wcs, W, H, i, unique, ra1,ra2,dec1,dec2)
            #print('Pixel', i, ': nu/ntot', nu, ntot)
            if nu > 0:
                i -= step
                break
            unique[:i,:] = False
            unique[H-1-i:,:] = False
            unique[:,:i] = False
            unique[:,W-1-i:] = False

        for j in range(max(i+1, 0), W//2):
            nu,ntot = self._ring_unique(wcs, W, H, j, unique, ra1,ra2,dec1,dec2)
            #print('Pixel', j, ': nu/ntot', nu, ntot)
            if nu == ntot:
                break
        return unique

    def _ring_unique(self,wcs, W, H, i, unique, ra1,ra2,dec1,dec2):
        lo, hix, hiy = i, W-i-1, H-i-1
        # one slice per side; we double-count the last pix of each side.
        sidex = slice(lo,hix+1)
        sidey = slice(lo,hiy+1)
        top = (lo, sidex)
        bot = (hiy, sidex)
        left  = (sidey, lo)
        right = (sidey, hix)
        xx = np.arange(W)
        yy = np.arange(H)
        nu,ntot = 0,0
        for slc in [top, bot, left, right]:
            #print('xx,yy', xx[slc], yy[slc])
            (yslc,xslc) = slc
            rr,dd = wcs.pixelxy2radec(xx[xslc]+1, yy[yslc]+1)
            U = (rr >= ra1 ) * (rr < ra2 ) * (dd >= dec1) * (dd < dec2)
            #print('Pixel', i, ':', np.sum(U), 'of', len(U), 'pixels are unique')
            unique[slc] = U
            nu += np.sum(U)
            ntot += len(U)
        #if allin:
        #    print('Scanned to pixel', i)
        #    break
        return nu,ntot

def main(doWhat=None,nproc=1,data_dir='./',dataset=None,date='mm-dd-yyyy',
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

    for cnt,brick in enumerate(bricks):
        if (cnt+1) % 10 == 0: 
            print('rank %d: %d/%d' % (comm.rank,cnt+1,len(bricks)))
        dr= derived_field_dir(brick,data_dir,date)
        try:
            os.makedirs(dr)
        except OSError:
            pass
        if doWhat == 'randoms_table':
            R= RandomsTables(data_dir,dataset,date=date)
            R.uniform_obiwana_obiwanb(brick)
        elif doWhat == 'heatmap_table':
            H= HeatmapTable(data_dir,dataset,date=date)
            H.write_table_for_datarelease(brick)
            

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

