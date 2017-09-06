from __future__ import print_function
import numpy as np
import fitsio
from glob import glob
import os
from astropy import units
from astropy.coordinates import SkyCoord
from astrometry.libkd.spherematch import match_radec

from astrometry.util.fits import fits_table, merge_tables
from tractor.brightness import NanoMaggies

def read_lines(fn_list):
    fin=open(fn_list,'r')
    lines=fin.readlines()
    fin.close()
    return np.sort(np.array( list(np.char.strip(lines)) ))

class ConvertTractor(object):
    """Convert DR3 fits_table to DR4 and beyond naming scheme
    """
    def isOld(self,cols):
        """cols: T.get_columns()"""
        if 'decam_flux' in cols:
            return True
        else:
            return False
   
    def old_to_new(self,T):
        """T: fits_table()
        """
        #T = fits_table(fn)
        #hdr = T.get_header()
        #primhdr = fitsio.read_header(fn)
        #os.environ['DUST_DIR']= '/home/kaylan/myrepo/dust'
        #outfn=fn.replace('.fits','-newformat.fits')
        return self.format_catalog(T, 'ugrizY', outfn=None,
                                   in_flux_prefix="decam_",
                                   flux_prefix='',write=False)
    
    def format_catalog(self, T, allbands, outfn,
                       in_flux_prefix='', flux_prefix='',write=False):
        # args: hdr, primhdr
        
        # Retrieve the bands in this catalog.
        bands = allbands
        print('Bands in this catalog:', bands)
        
        has_wise =    'wise_flux'    in T.columns()
        has_wise_lc = 'wise_lc_flux' in T.columns()
        has_ap =      'decam_apflux'       in T.columns()
        # Expand out FLUX and related fields from grz arrays to 'allbands'
        # (eg, ugrizY) arrays.
        B = np.array([allbands.index(band) for band in bands])
        keys = ['flux', 'flux_ivar', 'rchi2', 'fracflux', 'fracmasked', 'fracin',
                'nobs', 'anymask', 'allmask', 'psfsize', 'depth', 'galdepth']
        if has_ap:
            keys.extend(['apflux', 'apflux_resid', 'apflux_ivar'])    
            
        for k in keys:
            incol = '%s%s' % (in_flux_prefix, k)
            X = T.get(incol)
            # Handle array columns (eg, apflux)
            sh = X.shape
            if len(sh) == 3:
                nt,nb,N = sh
                A = np.zeros((len(T), len(allbands), N), X.dtype)
                A[:,B,:] = X
            else:
                A = np.zeros((len(T), len(allbands)), X.dtype)
                # If there is only one band, these can show up as scalar arrays.
                if len(sh) == 1:
                    A[:,B] = X[:,np.newaxis]
                else:
                    A[:,B] = X
            T.delete_column(incol)
            # FLUX_b for each band, rather than array columns.
            for i,b in enumerate(allbands):
                T.set('%s%s_%s' % (flux_prefix, k, b), A[:,i])
        
        # WISE
        wise_prefix='wise_'
        wise_keys=[]
        if has_wise:
            wise_keys.extend(['flux', 'flux_ivar','nobs','rchi2','fracflux'])
        if has_wise_lc:
            wise_keys.extend(['lc_flux','lc_flux_ivar','lc_nobs','lc_rchi2',
                              'lc_fracflux','lc_mjd'])
        for k in wise_keys:
            if 'lc_' in k:
                wisebands=['w1','w2']
            else:
                wisebands=['w1','w2','w3','w4']
            B_wise = np.array([wisebands.index(band) for band in wisebands])
            incol = '%s%s' % (wise_prefix, k)
            #raise ValueError
            X = T.get(incol)
            # Handle array columns (eg, apflux)
            sh = X.shape
            if len(sh) == 3:
                nt,nb,N = sh
                A = np.zeros((len(T), len(wisebands), N), X.dtype)
                A[:,B_wise,:] = X
            else:
                A = np.zeros((len(T), len(wisebands)), X.dtype)
                # If there is only one band, these can show up as scalar arrays.
                if len(sh) == 1:
                    A[:,B_wise] = X[:,np.newaxis]
                else:
                    A[:,B_wise]= X
            T.delete_column(incol)
            # FLUX_b for each band, rather than array columns.
            for i,b in enumerate(wisebands):
                T.set('%s%s_%s' % ('', k, b), A[:,i])
        ### Done WISE
        
        
        from tractor.sfd import SFDMap
        print('Reading SFD maps...')
        sfd = SFDMap()
        filts = ['%s %s' % ('DES', f) for f in allbands]
        wisebands = ['WISE W1', 'WISE W2', 'WISE W3', 'WISE W4']
        ebv,ext = sfd.extinction(filts + wisebands, T.ra, T.dec, get_ebv=True)
        T.ebv = ebv.astype(np.float32)
        ext = ext.astype(np.float32)
        decam_ext = ext[:,:len(allbands)]
        if has_wise:
            wise_ext  = ext[:,len(allbands):]
            
        wbands = ['w1','w2','w3','w4']
        
        trans_cols_opt  = []
        trans_cols_wise = []

        # No MW_TRANSMISSION_* columns at all
        for i,b in enumerate(allbands):
            col = 'mw_transmission_%s' % b
            T.set(col, 10.**(-decam_ext[:,i] / 2.5))
            trans_cols_opt.append(col)
            if has_wise:
                for i,b in enumerate(wbands):
                    col = 'mw_transmission_%s' % b
                    T.set(col, 10.**(-wise_ext[:,i] / 2.5))
                    trans_cols_wise.append(col)
                    
        from legacypipe.survey import release_number
        T.release = np.zeros(len(T), np.int16) + release_number
        
        # Column ordering...
        cols = ['release', 'brickid', 'brickname', 'objid', 'brick_primary', 
                'type', 'ra', 'dec', 'ra_ivar', 'dec_ivar',
                'bx', 'by', 'dchisq', 'ebv', 'mjd_min', 'mjd_max']
        
        def add_fluxlike(c):
            for b in allbands:
                cols.append('%s%s_%s' % (flux_prefix, c, b))
        def add_wiselike(c, bands=wbands):
            for b in bands:
                cols.append('%s_%s' % (c, b))
        
        add_fluxlike('flux')
        if has_wise:
            add_wiselike('flux')
            add_fluxlike('flux_ivar')
        if has_wise:
            add_wiselike('flux_ivar')
        if has_ap:
            for c in ['apflux', 'apflux_resid','apflux_ivar']:
                add_fluxlike(c)
                
        cols.extend(trans_cols_opt)
        cols.extend(trans_cols_wise)
        
        for c in ['nobs', 'rchisq', 'fracflux']:
            add_fluxlike(c)
            if has_wise:
                add_wiselike(c)
        for c in ['fracmasked', 'fracin', 'anymask', 'allmask']:
            add_fluxlike(c)
        if has_wise:
            pass # Only DR4 on has a wise mask
        #for i,b in enumerate(wbands[:2]):
        #    col = 'wisemask_%s' % (b)
        #    T.set(col, T.wise_mask[:,i])
        #    cols.append(col)
        for c in ['psfsize', 'psfdepth', 'galdepth']:
            add_fluxlike(c)
        
        if has_wise:
            cols.append('wise_coadd_id')
        if has_wise_lc:
            for c in ['lc_flux', 'lc_flux_ivar', 'lc_nobs', 'lc_fracflux',
                      'lc_rchisq','lc_mjd']:
                add_wiselike(c, bands=['w1','w2'])
                cols.extend([
                    'fracdev', 'fracdev_ivar',
                    'shapeexp_r', 'shapeexp_r_ivar',
                    'shapeexp_e1', 'shapeexp_e1_ivar',
                    'shapeexp_e2', 'shapeexp_e2_ivar',
                    'shapedev_r',  'shapedev_r_ivar',
                    'shapedev_e1', 'shapedev_e1_ivar',
                    'shapedev_e2', 'shapedev_e2_ivar',])
        # match case to T.
        cc = T.get_columns()
        cclower = [c.lower() for c in cc]
        for i,c in enumerate(cols):
            if (not c in cc) and c in cclower:
                j = cclower.index(c)
                cols[i] = cc[j]

        # Units
        deg = 'deg'
        degiv = '1/deg^2'
        arcsec = 'arcsec'
        flux = 'nanomaggy'
        fluxiv = '1/nanomaggy^2'
        units = dict(
            ra=deg, dec=deg, ra_ivar=degiv, dec_ivar=degiv, ebv='mag',
            shapeexp_r=arcsec, shapeexp_r_ivar='1/arcsec^2',
            shapedev_r=arcsec, shapedev_r_ivar='1/arcsec^2')
        # WISE fields
        wunits = dict(flux=flux, flux_ivar=fluxiv,
                      lc_flux=flux, lc_flux_ivar=fluxiv)
        # Fields that take prefixes (and have bands)
        funits = dict(
            flux=flux, flux_ivar=fluxiv,
            apflux=flux, apflux_ivar=fluxiv, apflux_resid=flux,
            psfdepth=fluxiv, galdepth=fluxiv, psfsize=arcsec)
        # add prefixes
        units.update([('%s%s' % (flux_prefix, k), v) for k,v in funits.items()])
        # add bands
        for b in allbands:
            units.update([('%s%s_%s' % (flux_prefix, k, b), v)
                          for k,v in funits.items()])
            # add WISE bands
        for b in wbands:
            units.update([('%s_%s' % (k, b), v)
                          for k,v in wunits.items()])
            
        # Create a list of units aligned with 'cols'
        units = [units.get(c, '') for c in cols]
        
        # Cleanup
        # Catalogue
        for k_new,k_old in zip(['psfdepth','rchisq'],
                               ['depth','rchi2']):
            for b in bands:
                T.rename('%s_%s' % (k_old,b),'%s_%s' % (k_new,b))
                #old_col= '%s_%s' % (k_old,b)
                #if old_col in T.get_columns():
                #    T.rename(old_col,'%s_%s' % (k_new,b))
                #else:
                #    print('warning'.upper() + ' %s not in fits_table' % old_col)
                    
        for k_new,k_old in zip(['rchisq'],
                               ['rchi2']):
            for b in ['w1','w2','w3','w4']:
                T.rename('%s_%s' % (k_old,b),'%s_%s' % (k_new,b))
        for k_new,k_old in zip(['lc_rchisq'],
                               ['lc_rchi2']):
            for b in ['w1','w2']:
                T.rename('%s_%s' % (k_old,b),'%s_%s' % (k_new,b))    
        for d_key in ['decam_mw_transmission']:
            if d_key in T.get_columns():
                T.delete_column(d_key)
        # ONLY in DR4 on
        for col in ['mjd_min','mjd_max','wise_coadd_id']:
            if col in cols:
                i= cols.index(col)
                cols.pop(i)
                units.pop(i)
        # check that cleanup worked
        for col in cols:
            if not col in T.get_columns():
                raise KeyError('%s in col but not in T' % col)
            
        #T.writeto(outfn, columns=cols, header=hdr, primheader=primhdr, units=units)
        if write:
            T.writeto(outfn, columns=cols, units=units)
            print('Wrote %s' % outfn)
        else:
            print('Skipping writing of %s' % outfn)
        return T



class CatalogueFuncs(object):
    '''funcs for using Dustins fits_table objects'''
    def stack(self,fn_list,textfile=True,
              shuffle=None):
        '''concatenates fits tables
        shuffle: set to an integer to randomly reads up to the first "shuffle" cats only
        '''
        if shuffle:
            assert( isinstance(shuffle, int))
        if textfile: 
            fns=read_lines(fn_list)
        else:
            fns= fn_list
        if len(fns) < 1: raise ValueError('Error: fns=',fns)
        if shuffle:
            print('shuffling %d' % shuffle)
            seed=7
            np.random.seed(seed)
            inds= np.arange(len(fns)) 
            np.random.shuffle(inds) 
            fns= fns[inds]
        cats= []
        for i,fn in enumerate(fns):
            print('reading %s %d/%d' % (fn,i+1,len(fns)))
            if shuffle and i >= shuffle: 
                print('shuffle_1000 turned ON, stopping read') 
                break 
            try:
                tab= fits_table(fn) 
                cats.append( tab )
            except IOError:
                print('Fits file does not exist: %s' % fn)
        return merge_tables(cats, columns='fillzero')

    def set_mags(self,cat):
        '''adds columns to fits_table cat'''
        # Remove white spaces
        cat.set('type', np.char.strip(cat.get('type')))
        # AB mags
        # Two kinds, as observed (including dust) and instrinsic (dust removed)
        for whichmag in ['wdust','nodust']:
            # DECam
            for band in ['g','r','z','w1','w2']:
                if whichmag == 'wdust':
                    flux= cat.get('flux_%s' % band)
                    flux_ivar= cat.get('flux_ivar_%s' % band)
                elif whichmag == 'nodust':
                    flux= cat.get('flux_%s' % band)/cat.get('mw_transmission_%s' % band)
                    flux_ivar= cat.get('flux_ivar_%s' % band)*\
                                np.power(cat.get('mw_transmission_%s' % band),2)
                else: raise ValueError()
                mag,mag_err= NanoMaggies.fluxErrorsToMagErrors(flux, flux_ivar)
                cat.set('mag_%s_%s' % (whichmag,band),mag)
                cat.set('mag_ivar_%s_%s' % (whichmag,band),1./np.power(mag_err,2))
        # Instrinsic fluxes
        whichmag='nodust'
        # DECam
        for band in ['g','r','z','w1','w2']:
            flux= cat.get('flux_%s' % band)/cat.get('mw_transmission_%s' % band)
            flux_ivar= cat.get('flux_ivar_%s' % band)*\
                                    np.power(cat.get('mw_transmission_%s' % band),2)
        cat.set('flux_%s_%s' % (whichmag,band),flux)
        cat.set('flux_ivar_%s_%s' % (whichmag,band),flux_ivar)

    def set_mags_OldDataModel(self,cat):
        '''for tractor catalogues with columns like decam_flux.shape(many,6)
        adds columns to fits_table cat'''
        # Remove white spaces
        cat.set('type', np.char.strip(cat.get('type')))
        # AB mags
        # Two kinds, as observed (including dust) and instrinsic (dust removed)
        for whichmag in ['wdust','nodust']:
            # DECam
            shp=cat.get('decam_flux').shape
            mag,mag_err= np.zeros(shp),np.zeros(shp)
            for iband in range(shp[1]):
                if whichmag == 'wdust':
                    flux= cat.get('decam_flux')[:,iband]
                    flux_ivar= cat.get('decam_flux_ivar')[:,iband]
                elif whichmag == 'nodust':
                    flux= cat.get('decam_flux')[:,iband]/cat.get('decam_mw_transmission')[:,iband]
                    flux_ivar= cat.get('decam_flux_ivar')[:,iband]*\
                                np.power(cat.get('decam_mw_transmission')[:,iband],2)
                else: raise ValueError()
                mag[:,iband],mag_err[:,iband]=NanoMaggies.fluxErrorsToMagErrors(flux, flux_ivar)
            cat.set('decam_mag_%s' % whichmag,mag)
            cat.set('decam_mag_ivar_%s' % whichmag,1./np.power(mag_err,2))
            # WISE
            if 'wise_flux' in cat.get_columns(): 
                shp=cat.get('wise_flux').shape
                mag,mag_err= np.zeros(shp),np.zeros(shp)
                for iband in range(shp[1]):
                    if whichmag == 'wdust':
                        flux= cat.get('wise_flux')[:,iband]
                        flux_ivar= cat.get('wise_flux_ivar')[:,iband]
                    elif whichmag == 'nodust':
                        flux= cat.get('wise_flux')[:,iband]/cat.get('wise_mw_transmission')[:,iband]
                        flux_ivar= cat.get('wise_flux_ivar')[:,iband]*\
                                    np.power(cat.get('wise_mw_transmission')[:,iband],2)
                    mag[:,iband],mag_err[:,iband]=NanoMaggies.fluxErrorsToMagErrors(flux, flux_ivar)
                cat.set('wise_mag_%s' % whichmag,mag)
                cat.set('wise_mag_ivar_%s' % whichmag,1./np.power(mag_err,2))
        # Instrinsic fluxes
        whichmag='nodust'
        # DECam
        shp=cat.get('decam_flux').shape
        flux,flux_ivar= np.zeros(shp),np.zeros(shp)
        for iband in range(shp[1]):
            flux[:,iband]= cat.get('decam_flux')[:,iband]/cat.get('decam_mw_transmission')[:,iband]
            flux_ivar[:,iband]= cat.get('decam_flux_ivar')[:,iband]*\
                                    np.power(cat.get('decam_mw_transmission')[:,iband],2)
        cat.set('decam_flux_%s' % whichmag,flux)
        cat.set('decam_flux_ivar_%s' % whichmag,flux_ivar)
        # WISE 
        if 'wise_flux' in cat.get_columns(): 
            shp=cat.get('wise_flux').shape
            flux,flux_err= np.zeros(shp),np.zeros(shp)
            for iband in range(shp[1]):
                flux[:,iband]= cat.get('wise_flux')[:,iband]/cat.get('wise_mw_transmission')[:,iband]
                flux_ivar[:,iband]= cat.get('wise_flux_ivar')[:,iband]*\
                                        np.power(cat.get('wise_mw_transmission')[:,iband],2)
            cat.set('wise_flux_%s' % whichmag,flux)
            cat.set('wise_flux_ivar_%s' % whichmag,flux_ivar)

        
class Cuts4MatchedCats(object):
    '''Cuts for MATCHED cats only'''
    def __init__(self,matched1,matched2):
        self.psf1 = (matched1.get('type') == 'PSF')
        self.psf2 = (matched2.get('type') == 'PSF')
        # Band dependent
        bands='grz'
        self.good={}
        for band,iband in zip(bands,range(6)):
            self.good[band]= ((matched1.get('flux_ivar_%s' % band) > 0) *\
                              (matched2.get('flux_ivar_%s' % band) > 0))      
        

#     def get_mags(self,cat):
#         mag= cat.get('decam_flux')/cat.get('decam_mw_transmission')
#         return 22.5 -2.5*np.log10(mag)

#     def get_mags_ivar(self,cat):
#         return np.power(np.log(10.)/2.5*cat.get('decam_flux'), 2)* \
#                                     cat.get('decam_flux_ivar')

class Matcher(object):
    '''sphere matches two ra,dec lists,
    ref,obs are astronometry.net "fits_tables" objects
    '''
    def __init__(self):
        pass
    def match_within(self,ref,obs,dist=1./3600):
        '''Find obs ra,dec that are within dist of ref ra,dec
        default=1 arcsec'''
        # Return 4 index arrays for indices where matches and where missing
        imatch=dict(ref=[],obs=[])
        imiss=dict(ref=[],obs=[])
        # cat1 --> ref cat
        # cat2 --> cat matching to the ref cat
        #if False:
        #    cat1 = SkyCoord(ra=ref.get('ra')*units.degree, dec=ref.get('dec')*units.degree)
        #    cat2 = SkyCoord(ra=obs.get('ra')*units.degree, dec=obs.get('dec')*units.degree)
        #    idx, d2d, d3d = cat1.match_to_catalog_3d(cat2)
        #    b= np.array(d2d) <= dist
        #    imatch['ref']= np.arange(len(ref))[b]
        #    imatch['obs']= np.array(idx)[b]
        I,J,d2d = match_radec(ref.ra,ref.dec, obs.ra,obs.dec, dist,
                              nearest=True)
        imatch['ref']= I
        imatch['obs']= J
        # 
        print("Matched: %d/%d objects" % (imatch['ref'].size,len(ref)))
        imiss['ref'] = np.delete(np.arange(len(ref)), imatch['ref'], axis=0)
        imiss['obs'] = np.delete(np.arange(len(obs)), imatch['obs'], axis=0)
        return imatch,imiss,d2d

    def nearest_neighbor(self,ref):
        '''return indices of nearest nearsest and distances to them'''
        cat1 = SkyCoord(ra=ref.get('ra')*units.degree, dec=ref.get('dec')*units.degree)
        cat2 = SkyCoord(ra=ref.get('ra')*units.degree, dec=ref.get('dec')*units.degree)
        idx, d2d, d3d = cat1.match_to_catalog_3d(cat2,nthneighbor=2)
        #b= np.array(d2d) <= within
        ref_nn= np.arange(len(ref))
        obs_nn= np.array(idx)
        dist= np.array(d2d)
        return ref_nn,obs_nn,dist

    def nearest_neighbors_within(self,ref,obs,within=1./3600,min_nn=1,max_nn=5):
        '''Find obs ra,dec that are within dist of ref ra,dec
        default=1 arcsec'''
        # cat1 --> ref cat
        # cat2 --> cat matching to the ref cat
        cat1 = SkyCoord(ra=ref.get('ra')*units.degree, dec=ref.get('dec')*units.degree)
        cat2 = SkyCoord(ra=obs.get('ra')*units.degree, dec=obs.get('dec')*units.degree)
        ref_nn,obs_nn,dist={},{},{}
        for nn in range(min_nn,max_nn+1):
            idx, d2d, d3d = cat1.match_to_catalog_3d(cat2,nthneighbor=nn)
            b= np.array(d2d) <= within
            ref_nn[str(nn)]= np.arange(len(ref))[b]
            obs_nn[str(nn)]= np.array(idx)[b]
            dist[str(nn)]= np.array(d2d)[b]
            print("within 1arcsec, nn=%d, %d/%d" % (nn,ref_nn[str(nn)].size,len(ref)))
        return ref_nn,obs_nn,dist

class TargetTruth(object):
    '''Build Target Truth catalogues, matching to DR3
    '''
    def __init__(self):
        self.truth_dir='/project/projectdirs/desi/target/analysis/truth'
        self.dr3_dir='/global/project/projectdirs/cosmo/data/legacysurvey/dr3'
        self.save_dir='/project/projectdirs/desi/users/burleigh/desi/target/analysis/truth'
  
    def in_region(self,ra,dec, rlo=0.,rhi=360.,dlo=0.,dhi=30.):
        '''ra,dec are numpy arrays,
        return indices of ra,dec where they are inside specified region'''
        if rlo < rhi:
            return (ra >= rlo) * (ra <= rhi) *\
                   (dec >= dlo) * (dec <= dhi)
        else: # RA wrap
            return np.logical_or(ra >= rlo, ra <= rhi) *\
                   (dec >= dlo) * (dec <= dhi)

    def bricks_in_region(self,rlo=0.,rhi=360.,dlo=0.,dhi=30.):
        print('Region: rlo=%.2f, rhi=%.2f, dlo=%.2f, dhi=%.2f' % (rlo,rhi,dlo,dhi))
        bricks=fits_table(os.path.join(self.dr3_dir,'survey-bricks.fits.gz'))
        i={}
        # Loop over 4 corners of each Brick
        for cnt,(ra,dec) in zip(range(1,5),[('ra1','dec1'),('ra1','dec2'),('ra2','dec1'),('ra2','dec2')]):
            i[str(cnt)]= self.in_region(bricks.get(ra),bricks.get(dec),\
                                        rlo=rlo,rhi=rhi,dlo=dlo,dhi=dhi)
            print('corner=%s, number bricks=%d' % (str(cnt),len(bricks.get('ra')[ i[str(cnt)] ])))
        i= np.any((i['1'],i['2'],i['3'],i['4']),axis=0)
        print('any corner, number bricks=%d' % len(bricks.get('ra')[i]))
        names= bricks.get('brickname')[i] 
        if not len(list(set(names))) == len(names):
            raise ValueError('Repeated brick names')
        return names

    def sweep_corner2radec(self,text='000m005'):
        ra,dec= float(text[:3]),float(text[4:])
        if text[3] == 'm': 
            dec*= -1
        return ra,dec

    def sweeps_in_region(self,rlo=0.,rhi=360.,dlo=0.,dhi=30.):
        print('Region: rlo=%.2f, rhi=%.2f, dlo=%.2f, dhi=%.2f' % (rlo,rhi,dlo,dhi))
        fns=glob(os.path.join(self.dr3_dir,'sweep/3.0/','sweep-*.fits'))
        fns=np.array(fns)
        assert(len(fns) > 0)
        # Loop over sweep regions
        # Ask if any corner of each region is in the region we are interested in
        keep=np.zeros(len(fns)).astype(bool)
        for cnt,fn in enumerate(fns):
            left,right= os.path.basename(fn).split('.')[0].split('-')[1:]
            ra1,dec1= self.sweep_corner2radec(text=left)
            ra2,dec2= self.sweep_corner2radec(text=right)
            # Loop over the 4 corners
            b=False
            for ra,dec in [(ra1,dec1),(ra1,dec2),(ra2,dec1),(ra2,dec2)]:
                b= np.logical_or(b, self.in_region(ra,dec,\
                                                   rlo=rlo,rhi=rhi,dlo=dlo,dhi=dhi)
                                 )
            if b:
                print(ra1,'-->',ra2,'  ',dec1,'-->',dec2)
                keep[cnt]= True
        if not len(fns[keep]) > 0:
            raise ValueError('Something amiss, no sweeps overlap region')
        print('%d/%d sweeps in the region' % (len(fns[keep]),len(fns)))
        return fns[keep] 

    def cosmos_zphot(self):
        # Data
        cosmos=fits_table(os.path.join(self.truth_dir,'cosmos-zphot.fits.gz'))
        # Bricks
        bnames= self.bricks_in_region(rlo=cosmos.get('ra').min(), rhi=cosmos.get('ra').max(),\
                                          dlo=cosmos.get('dec').min(),dhi=cosmos.get('dec').max())
        # Tractor Catalogues --> file list
        catlist= os.path.join(self.save_dir,'cosmos_dr3_bricks.txt')
        if not os.path.exists(catlist):
            fout=open(catlist,'w')
            for b in bnames:
                fn= os.path.join(self.dr3_dir,'tractor/%s/tractor-%s.fits' % (b[:3],b))
                fout.write('%s\n' % fn)
            fout.close()
            print('Wrote %s' % catlist)
        # Match
        fits_funcs= CatalogueFuncs()
        dr3=fits_funcs.stack(os.path.join(self.save_dir,'cosmos_dr3_bricks.txt'))
        mat=Matcher()
        imatch,imiss,d2d= mat.match_within(cosmos,dr3) #,dist=1./3600)
        cosmos.cut(imatch['ref'])
        dr3.cut(imatch['obs'])
        # Save
        cosmos.writeto(os.path.join(self.save_dir,'cosmos-zphot-dr3matched.fits'))
        dr3.writeto(os.path.join(self.save_dir,'dr3-cosmoszphotmatched.fits'))
        print('Wrote %s\nWrote %s' % (os.path.join(self.save_dir,'cosmos-zphot-dr3matched.fits'),\
                                      os.path.join(self.save_dir,'dr3-cosmoszphotmatched.fits')))

    def vipers(self):
        # Data
        w1=fits_table(os.path.join(self.truth_dir,'vipers-w1.fits.gz'))
        w4=fits_table(os.path.join(self.truth_dir,'vipers-w4.fits.gz'))
        # Bricks
        for data in [w1,w4]:
            data.set('ra',data.get('alpha'))
            data.set('dec',data.get('delta'))
        bnames={}
        for data,key in zip([w1,w4],['w1','w4']):
            bnames[key]= self.bricks_in_region(rlo=data.get('ra').min(), rhi=data.get('ra').max(),\
                                              dlo=data.get('dec').min(),dhi=data.get('dec').max())
        bricks=np.array([])
        for key in bnames.keys():
            bricks=np.concatenate((bricks,bnames[key]))
        # Tractor Catalogues --> file list
        catlist= os.path.join(self.save_dir,'vipers_dr3_bricks.txt')
        if not os.path.exists(catlist):
            fout=open(catlist,'w')
            for b in bricks:
                fn= os.path.join(self.dr3_dir,'tractor/%s/tractor-%s.fits' % (b[:3],b))
                fout.write('%s\n' % fn)
            fout.close()
            print('Wrote %s' % catlist)
        # Merge w1,w4 for matching
        vipers= []
        for fn in [os.path.join(self.truth_dir,'vipers-w1.fits.gz'),\
                   os.path.join(self.truth_dir,'vipers-w4.fits.gz')]:
            vipers.append( fits_table(fn) )
        vipers= merge_tables(vipers, columns='fillzero')
        vipers.set('ra',data.get('alpha'))
        vipers.set('dec',data.get('delta'))
        # Match
        fits_funcs= CatalogueFuncs()
        dr3=fits_funcs.stack(os.path.join(self.save_dir,'vipers_dr3_bricks.txt'))
        mat=Matcher()
        imatch,imiss,d2d= mat.match_within(vipers,dr3) #,dist=1./3600)
        vipers.cut(imatch['ref'])
        dr3.cut(imatch['obs'])
        # Save
        vipers.writeto(os.path.join(self.save_dir,'vipersw1w4-dr3matched.fits'))
        dr3.writeto(os.path.join(self.save_dir,'dr3-vipersw1w4matched.fits'))
        print('Wrote %s\nWrote %s' % (os.path.join(self.save_dir,'vipersw1w4-dr3matched.fits'),\
                                      os.path.join(self.save_dir,'dr3-vipersw1w4matched.fits')))

    def deep2(self):
        # Data
        deep2={}
        for key in ['1','2','3','4']:
            deep2[key]=fits_table('/project/projectdirs/desi/target/analysis/truth/deep2-field%s.fits.gz' % key)
        # Bricks
        bnames={}
        for key in deep2.keys():
            bnames[key]= self.bricks_in_region(rlo=deep2[key].get('ra').min(), rhi=deep2[key].get('ra').max(),\
                                              dlo=deep2[key].get('dec').min(),dhi=deep2[key].get('dec').max())
            print('Field=%s, Num Bricks=%d, Bricks:' % (key,len(bnames[key])), bnames[key])
        bricks=np.array([])
        for key in bnames.keys():
            bricks=np.concatenate((bricks,bnames[key]))
        # Tractor Catalogues --> file list
        catlist= os.path.join(self.save_dir,'deep2_dr3_bricks.txt')
        if not os.path.exists(catlist):
            fout=open(catlist,'w')
            for b in bricks:
                fn= os.path.join(self.dr3_dir,'tractor/%s/tractor-%s.fits' % (b[:3],b))
                fout.write('%s\n' % fn)
            fout.close()
            print('Wrote %s' % catlist)
        # Merge for matching
        dp2= [deep2['2'],deep2['3'],deep2['4']]
        dp2= merge_tables(dp2, columns='fillzero')
        # Match
        fits_funcs= CatalogueFuncs()
        dr3=fits_funcs.stack(os.path.join(self.save_dir,'deep2_dr3_bricks.txt'))
        mat=Matcher()
        imatch,imiss,d2d= mat.match_within(dp2,dr3) #,dist=1./3600)
        dp2.cut(imatch['ref'])
        dr3.cut(imatch['obs'])
        fits_funcs.set_extra_data(dr3)
        # Save
        dp2.writeto(os.path.join(self.save_dir,'deep2f234-dr3matched.fits'))
        dr3.writeto(os.path.join(self.save_dir,'dr3-deep2f234matched.fits'))
        print('Wrote %s\nWrote %s' % (os.path.join(self.save_dir,'deep2f234-dr3matched.fits'),\
                                      os.path.join(self.save_dir,'dr3-deep2f234matched.fits')))

    def qso(self):
        # Christophe's qso catalogue
        qso=fits_table(os.path.join(self.save_dir,'CatalogQSO.fits.gz'))
        # DR7 and boss only
        keep={}
        for key in list(set(qso.get('source'))):
            print('%s: %d' % (key,len(qso[qso.get('source') == key])))
            keep[key.strip()]= qso.get('source') == key
        qso.cut( np.logical_or(keep['BOSS'],keep['SDSSDR7QSO']) )
        # Stripe 82
        rlo,rhi= 315., 45.
        dlo,dhi= -1.25, 1.25
        qso.cut( self.in_region(qso.get('ra'),qso.get('dec'), \
                                rlo=rlo,rhi=rhi,dlo=dlo,dhi=dhi) )
        # Bricks
        sweeps= self.sweeps_in_region(rlo=rlo, rhi=rhi,\
                                      dlo=dlo,dhi=dhi)
        print('sweeps= ',sweeps)
        sys.exit('early')
        # Tractor Catalogues --> file list
        #catlist= os.path.join(self.save_dir,'CatalogQSO_dr3_sweeps.txt')
        #if not os.path.exists(catlist):
        #    fout=open(catlist,'w')
        #    for b in bricks:
        #        fn= os.path.join(self.dr3_dir,'tractor/%s/tractor-%s.fits' % (b[:3],b))
        #        fout.write('%s\n' % fn)
        #    fout.close()
        #    print('Wrote %s' % catlist)
        # Match
        fits_funcs= CatalogueFuncs()
        #dr3=fits_funcs.stack(os.path.join(self.save_dir,'CatalogQSO_dr3_bricks.txt'))
        dr3=fits_funcs.stack(sweeps,textfile=False)
        mat=Matcher()
        imatch,imiss,d2d= mat.match_within(qso,dr3) #,dist=1./3600)
        qso.cut(imatch['ref'])
        dr3.cut(imatch['obs'])
        fits_funcs.set_extra_data(dr3)
        # Save
        qso.writeto(os.path.join(self.save_dir,'qso-dr3sweepmatched.fits'))
        dr3.writeto(os.path.join(self.save_dir,'dr3-qsosweepmatched.fits'))
        print('Wrote %s\nWrote %s' % (os.path.join(self.save_dir,'qso-dr3sweepmatched.fits'),\
                                      os.path.join(self.save_dir,'dr3-qsosweepmatched.fits')))

    def get_acs_six_col(self):
        savenm= os.path.join(self.save_dir,'acs_six_cols.fits')
        if not os.path.exists(savenm):
            acsfn=os.path.join(self.save_dir,'ACS-GC_published_catalogs',
                               'acs_public_galfit_catalog_V1.0.fits.gz')
            acs=fits_table(acsfn) 
            # Repackage
            tab= fits_table()
            for key in ['ra','dec','re_galfit_hi','n_galfit_hi',
                        'ba_galfit_hi','pa_galfit_hi']:
                tab.set(key, acs.get(key))
            # Cuts & clean
            tab.cut( tab.flag_galfit_hi == 0 )
            # -90,+90 --> 0,180 
            tab.pa_galfit_hi += 90. 
            # Save
            tab.writeto(savenm)
            print('Wrote %s' % savenm)
            # Save
            qso.writeto(os.path.join(self.save_dir,'qso-dr3sweepmatched.fits'))
            dr3.writeto(os.path.join(self.save_dir,'dr3-qsosweepmatched.fits'))
            print('Wrote %s\nWrote %s' % (os.path.join(self.save_dir,'qso-dr3sweepmatched.fits'),\
                                          os.path.join(self.save_dir,'dr3-qsosweepmatched.fits')))
