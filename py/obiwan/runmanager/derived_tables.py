"""
Using obiwan outputs, make a table of 'official' randoms per-brick
- uniform randoms: random ra,dec + geometry cut
- obiwan randoms: uniform randoms + recovered by tractor
"""

import numpy as np
import os
from glob import glob
import pandas as pd
from collections import Counter

from obiwan.db_tools import all_psqlcols_for_ids
try: 
    from astrometry.util.fits import fits_table, merge_tables
    from astrometry.util.util import Tan
    from astrometry.libkd.spherematch import match_radec
except ImportError:
    pass

def derived_field_dir(brick,data_dir,date):
    return os.path.join(data_dir,'derived_%s' % date,
                        brick[:3],brick)

def datarelease_dir(eboss_or_desi):
    proj='/global/project/projectdirs/cosmo/data/legacysurvey'
    if eboss_or_desi == 'eboss':
        dr= 'dr3'
    elif eboss_or_desi == 'desi':
        dr='dr5'
    return os.path.join(proj,dr)
       
def is_numeric(obj):                                        
    try: 
        tmp=obj+5
    except TypeError:
        return False
    return True

class Bit(object):
    def set(self,value, bit):
        """change bit to 1, bit is 0-indexed"""
        return value | (1<<bit)

    def clear(self,value, bit):
        """change bit to 0, bit is 0-indexed"""
        return value & ~(1<<bit)

class RandomsTable(object):
    """Creates the uniform,obiwan_a,obiwan_b randoms tables for a single brick

    Final table has same number rows as uniform table, and the obiwan
        rows are filled in wherever there is a matching unique_id
        between uniform and obiwan. A bitmask column 'obiwan_mask' which says
        whether the random was recovered by Tractor or not, and whether
        the random is near a previously existing real source in a DR 
        catalogue, like DR3 or DR5
    """
    def __init__(self, data_dir,eboss_or_desi,db_randoms_table,
                 date='mm-dd-yyyy'):
        self.data_dir= data_dir
        self.eboss_or_desi= eboss_or_desi
        self.db_randoms_table= db_randoms_table
        self.date= date

    def run(self,brick):
        tab= self.merge_randoms_tables(brick)
        self.add_flag_for_realsources(tab,brick)
        # Write
        derived_dir= derived_field_dir(brick,self.data_dir,self.date)
        fn= os.path.join(derived_dir,'randoms.fits')
        self.write_table(tab,fn)
     
    def merge_randoms_tables(self,brick):
        """Computes final joined randoms tables

        Includes uniform randoms, info from psql db, which of these were recovered
            by Tractor, and the associated tractor info for those 

        Args:
            brick: brickname
        
        Returns: 
            joined randoms table
        """
        search= os.path.join(self.data_dir,'tractor',
                             brick[:3],brick,
                             'rs*','tractor-%s.fits' % brick)
        rsdirs= glob(search)
        rsdirs= [os.path.dirname(dr)
                 for dr in rsdirs]
        if len(rsdirs) == 0:
            raise ValueError('no rsdirs found: %s' % search)
        uniform=[]
        for dr in rsdirs:
            simcat= fits_table((os.path.join(dr,'simcat-elg-%s.fits' % brick)
                                .replace('/tractor/','/obiwan/')))
            idsadded= fits_table((os.path.join(dr,'sim_ids_added.fits')
                                .replace('/tractor/','/obiwan/')))
            # Uniform randoms (injected at touching at least 1 ccd)
            assert(len(idsadded) == len(set(idsadded.id)))
            simcat.cut( pd.Series(simcat.id).isin(idsadded.id) )
            simcat.set('unique_id',self.unique_id(simcat.id.astype(str),
                                                  brick,os.path.basename(dr)))
            self.add_psql_to_uniform_table(simcat,self.db_randoms_table)
            # Recovered by Tractor
            tractor= fits_table(os.path.join(dr,'tractor-%s.fits' % brick))
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
            add_vals={}
            for trac_key in tractor.get_columns():
                key= 'tractor_'+trac_key
                if is_numeric(tractor.get(trac_key)):
                    shp= (len(simcat),) + tractor.get(trac_key).shape[1:]
                    add_vals[key]= np.zeros(shp) +np.nan
                else:
                    add_vals[key]= np.array(['']*len(simcat))
                add_vals[key][I]= tractor.get(trac_key)
                simcat.set(key,add_vals[key])
            # Mask
            mask= np.zeros(len(simcat),dtype=np.int8)
            mask[I]= Bit().set(mask[I],0)
            simcat.set('obiwan_mask',mask)
            # add to list uniform tables
            uniform.append(simcat)
        return merge_tables(uniform, columns='fillzero')

    def unique_id(self,id_array,brick,rs_dir):
        """For a given random injected into a brick during a given iteration

        Args:
            id_array: randoms ids 
            brick: brick
            rs_dir: like rs0 or rs300
        """
        ids= np.array(id_array,dtype=object) + "_%s_%s" % (brick,rs_dir)
        # FITS can't handle numpy type 'object'
        return ids.astype(str)

    def add_psql_to_uniform_table(self,uniform,db_randoms_table):
        """Add randoms db columns from psql to the uniform randoms table

        Args:
            uniform: fits table
            db_randoms_table: name of the psql db table
        """
        db_dict= all_psqlcols_for_ids(uniform.id, db_randoms_table=db_randoms_table)
        assert(all(db_dict['id'] - uniform.id == 0))
        for key,val in db_dict.items():
            if key in ['id']:
                pass
            uniform.set('psql_%s' % key,val)

    def add_flag_for_realsources(self,tab,brick):
        """Flag sources also in DR3, DR5

        Args:
            tab: table returned by merged_randoms_table()
        """
        real= fits_table(os.path.join(datarelease_dir(self.eboss_or_desi),
                                      'tractor',brick[:3],
                                      'tractor-%s.fits' % brick))
        # nearest match in (ra2,dec2) for each point in (ra1,dec1)
        I,J,d = match_radec(tab.ra,tab.dec,
                            real.ra,real.dec, 1./3600,
                            nearest=True)
        assert(np.all(d <= 1./3600))
        bool_matched= np.zeros(len(tab),bool)
        bool_matched[I]= True
        recovered_and_matched= ((tab.obiwan_mask == 1) & 
                                (bool_matched))
        if len(tab[recovered_and_matched]) > 0:
            mask= tab.obiwan_mask
            mask[recovered_and_matched]= Bit().set(mask[recovered_and_matched],1)
            tab.set('obiwan_mask',mask)

    def write_table(self,tab,fn):
        """Write the merged randoms table is doesn't already exist"""
        if not os.path.exists(fn): 
            tab.writeto(fn)
            print('Wrote %s' % fn)
    

class TargetSelection(object):
    def __init__(self,eboss_or_desi): 
        assert(eboss_or_desi in ['eboss','desi'])
        self.eboss_or_desi= eboss_or_desi

    def keep(self,tractor):
        if self.eboss_or_desi == 'eboss':
            return self.ebossIsElg(tractor)
        elif self.eboss_or_desi == 'desi':
            return self.desiIsElg(tractor)

    def desiIsElg(self,tractor):
        kw={}
        if 'brick_primary' in tractor.get_columns():
            kw.update(primary=tractor.brick_primary)
        for band,iband in [('g',1),('r',2),('z',4)]:
            kw[band+'flux']= tractor.get('flux_'+band) / tractor.get('mw_transmission_'+band)
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
        inRegion=dict(ngc= ((tractor.ra > 126.) &
                            (tractor.ra < 168.) &
                            (tractor.dec > 14.) &
                            (tractor.ra < 34.)),
                      sgc_a= ((tractor.ra > 317.) &
                              (tractor.ra < 360.) &
                              (tractor.dec > -2.) &
                              (tractor.ra < 2.)),
                      sgc_b= ((tractor.ra > 0.) &
                              (tractor.ra < 45.) &
                              (tractor.dec > -5.) &
                              (tractor.ra < 5.)))
        inRegion.update(sgc= ((inRegion['sgc_a']) | 
                              (inRegion['sgc_b'])))
        # tycho2inblob == False
        # SDSS bright object mask & 0 < V < 11.5 mag Tycho2 stars mask
        # anymask[grz] == 0
        # custom mask for eboss23
        self.add_grz_mag(tractor)
        gr= tractor.gmag - tractor.rmag
        rz= tractor.rmag - tractor.zmag
        colorCut= dict(sgc= ((tractor.gmag > 21.825) &
                             (tractor.gmag < 22.825) &
                             (-0.068 * rz + 0.457 < gr) &
                             (gr < 0.112 * rz + 0.773) &
                             (0.218 * gr + 0.571 < rz) &
                             (rz < -0.555 * gr + 1.901)),
                       ngc= ((tractor.gmag > 21.825) &
                             (tractor.gmag < 22.9) &
                             (-0.068 * rz + 0.457 < gr) &
                             (gr < 0.112 * rz + 0.773) &
                             (0.637 * gr + 0.399 < rz) &
                             (rz < -0.555 * gr + 1.901)))
        return ((inRegion['ngc'] & colorCut['ngc']) |
                (inRegion['sgc'] & colorCut['sgc']))

    def add_grz_mag(self,tractor):
        for band,iband in [('g',1),('r',2),('z',4)]:
            flux_ext= tractor.get('flux_'+band) / tractor.get('mw_transmission_'+band)
            tractor.set(band+'mag',self.flux2mag(flux_ext))

    def flux2mag(self,nmgy):
        return -2.5 * (np.log10(nmgy) - 9)


class TargetsTable(object):
    """Apply target selection to the RandomsTables and write it out per-brick"""
    def __init__(self,data_dir,eboss_or_desi,date='mm-dd-yyyy'):
        self.data_dir= data_dir
        self.eboss_or_desi= eboss_or_desi
        self.date= date
    
    def run(self,brick):
        derived_dir= derived_field_dir(brick,self.data_dir,self.date)
        for randoms_tab in ['uniform','obiwan_a','obiwan_b','obiwan_real']:
            self.write_targets(derived_dir,
                               randoms_table=randoms_tab)

    def write_targets(self,derived_dir,
                      randoms_table=None):
        TS= TargetSelection(self.eboss_or_desi) 
        fn= os.path.join(derived_dir,'randoms_%s.fits' % randoms_table)
        if randoms_table in ['uniform']:
            tractor= fits_table(fn)
            if 'gflux' in tractor.get_columns():
                for band in 'grz':
                    tractor.rename(band+'flux','flux_'+band)
        else:
            tractor= self.read_tractor(fn)
        tractor.cut(TS.keep(tractor))
        savefn= fn.replace('.fits','_%s.fits' % self.eboss_or_desi)
        tractor.writeto(savefn)
        print('Wrote %s' % savefn)

    def read_tractor(self,tractor_fn):
        columns=['brick_primary', 'type','ra','dec',
                 'brickname']
        for band in 'grz':
            for prefix in ['flux_','mw_transmission_',
                           'allmask_','anymask_']:
                columns.append(prefix+band)
        return fits_table(tractor_fn,columns=columns)




class HeatmapTable(object):
    """Create a fits table with the heatmap values, one per brick"""
    def __init__(self, data_dir,eboss_or_desi,date='mm-dd-yyyy'):
        self.data_dir= data_dir
        self.eboss_or_desi= eboss_or_desi
        self.date= date
        self.surveyBricks = fits_table(os.path.join(os.environ['LEGACY_SURVEY_DIR'],
                                                    'survey-bricks.fits.gz'))

    def run(self,brick):
        self.write_table_using_datarelease(brick)
        self.write_table_using_obiwan(brick)

    def write_table_using_datarelease(self,brick):
        derived_dir= derived_field_dir(brick,self.data_dir,self.date)
        fn= os.path.join(derived_dir,'heatmap_datarelease.fits')
        if os.path.exists(fn): 
            print('Skipping, already exist: ',fn)
        else:
            tab= self.get_table_for_datarelease([brick])
            if tab:
                tab.writeto(fn)
                print('Wrote %s' % fn)

    def write_table_using_obiwan(self,brick):
        derived_dir= derived_field_dir(brick,self.data_dir,self.date)
        fn= os.path.join(derived_dir,'heatmap_datarelease.fits')
        if os.path.exists(fn): 
            print('Skipping, already exist: ',fn)
        else:
            tab= self.get_table_for_datarelease([brick])
            if tab:
                tab.writeto(fn)
                print('Wrote %s' % fn)
 

    def get_table_for_datarelease(self,bricklist):
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
       
        dirprefix= datarelease_dir(self.eboss_or_desi)
        for ibrick,brick in enumerate(bricklist):
            #words = fn.split('/')
            #dirprefix = '/'.join(words[:-4])
            #print('Directory prefix:', dirprefix)
            #words = words[-4:]
            #brick = words[2]
            #print('Brick', brick)
            tfn = os.path.join(dirprefix, 'tractor', brick[:3], 'tractor-%s.fits'%brick)
            if self.eboss_or_desi == 'desi': # DR5 version of tractor cats
                columns=['brick_primary', 'type',
                         'psfsize_g', 'psfsize_r', 'psfsize_z',
                         'psfdepth_g', 'psfdepth_r', 'psfdepth_z',
                         'galdepth_g', 'galdepth_r', 'galdepth_z',
                         'ebv',
                         'mw_transmission_g', 'mw_transmission_r', 'mw_transmission_z',
                         'nobs_w1', 'nobs_w2', 'nobs_w3', 'nobs_w4',
                         'nobs_g', 'nobs_r', 'nobs_z',
                         'mw_transmission_w1', 'mw_transmission_w2', 'mw_transmission_w3', 'mw_transmission_w4']
            elif self.eboss_or_desi == 'eboss': # DR3 version of tractor cats
                columns=['brick_primary', 'type', 'decam_psfsize',
                         'decam_depth', 'decam_galdepth',
                         'ebv', 'decam_mw_transmission',
                         'decam_nobs',
                         'wise_nobs', 'wise_mw_transmission']
            try:
                T = fits_table(tfn, columns=columns)
                #print('Read %s' % tfn)
            except:
                print('Failed to read %s' % tfn)
                return None
                #print('Failed to read FITS table', tfn)
                #import traceback
                #traceback.print_exc()
                #print('Carrying on.')
                #continue


            if self.eboss_or_desi == 'desi':
                hasBands= [band for band in 'grz' if any(T.get('nobs_'+band) > 0)]
            elif self.eboss_or_desi == 'eboss':
                hasBands= [band 
                           for band,iband in [('g',1),('r',2),('z',4)]
                           if any(T.decam_nobs[:,iband] > 0)]
            brickset.add(brick)
            gn.append(0)
            rn.append(0)
            zn.append(0)

            gnhist.append([0 for i in range(nnhist)])
            rnhist.append([0 for i in range(nnhist)])
            znhist.append([0 for i in range(nnhist)])

            index = -1
            ibrick = np.nonzero(self.surveyBricks.brickname == brick)[0][0]
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

            if self.eboss_or_desi == 'desi':
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
                
            elif self.eboss_or_desi == 'eboss':
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

            br = self.surveyBricks[ibrick]

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
        
      
            # Does a check on the legacysurvey-{brick}-nexp*.fits files
            if False: 
                #filepart = words[-1]
                #filepart = filepart.replace('.fits.gz', '')
                #filepart = filepart.replace('.fits.fz', '')
                #print('File:', filepart)
                #band = filepart[-1]
                #assert(band in 'grz')
                nlist,nhist = dict(g=(gn,gnhist), r=(rn,rnhist), z=(zn,znhist))[band]
                for band in hasBands:
                    fn= os.path.join(dirprefix, 'coadd', 
                                     brick[:3],brick,
                                     'legacysurvey-%s-nexp-%s.fits.gz' % (brick,band))
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
        T.ra  = self.surveyBricks.ra [ibricks]
        T.dec = self.surveyBricks.dec[ibricks]
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

def main_mpi(bricks=[],doWhat=None,eboss_or_desi=None,
             db_randoms_table=None,
             nproc=1,data_dir='./',date='mm-dd-yyyy'):
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
        comm= MyComm()

    if doWhat == 'randoms':
        tabMaker= RandomsTable(data_dir,eboss_or_desi,db_randoms_table,
                               date=date)
    elif doWhat == 'targets':
        tabMaker= TargetsTable(data_dir,eboss_or_desi,date=date)
    elif doWhat == 'heatmap':
        tabMaker= HeatmapTable(data_dir,eboss_or_desi,date=date)
    
    for cnt,brick in enumerate(bricks):
        if (cnt+1) % 10 == 0: 
            print('rank %d: %d/%d' % (comm.rank,cnt+1,len(bricks)))
        dr= derived_field_dir(brick,data_dir,date)
        try:
            os.makedirs(dr)
        except OSError:
            pass
        tabMaker.run(brick)
            

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--doWhat', type=str, choices=['randoms','targets','heatmap'],required=True)
    parser.add_argument('--data_dir', type=str, required=True, 
                        help='path to obiwan/, tractor/ dirs') 
    parser.add_argument('--db_randoms_table', type=str, choices=['obiwan_eboss_elg',
                                    'obiwan_elg_dr5','obiwan_cosmos'],required=True)
    parser.add_argument('--nproc', type=int, default=1, help='set to > 1 to run mpi4py') 
    parser.add_argument('--bricks_fn', type=str, default=None,
                        help='specify a fn listing bricks to run, or a single default brick will be ran') 
    parser.add_argument('--eboss_or_desi', type=str, choices=['eboss','desi'], 
                        help='for obiwan_randoms_b',required=True) 
    parser.add_argument('--date', type=str,help='mm-dd-yyyy, to label derived directory by',required=True) 
    args = parser.parse_args()
    
    # Bricks to run
    if args.bricks_fn is None:
        bricks= ['1266p292']
    else:
        bricks= np.loadtxt(args.bricks_fn,dtype=str)

    kwargs= vars(args)
    for dropCol in ['bricks_fn']:
        del kwargs[dropCol]
    kwargs.update(bricks=bricks)

    main_mpi(**kwargs)

