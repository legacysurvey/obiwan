"""
Using obiwan outputs, make a table of 'official' randoms per-brick
- uniform randoms: random ra,dec + geometry cut
- obiwan randoms: uniform randoms + recovered by tractor
"""

import numpy as np
import os
from glob import glob
import pandas as pd
from collections import Counter,defaultdict

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

def datarelease_dir(drNumber):
    assert(drNumber in ['dr3','dr5'])
    proj='/global/project/projectdirs/cosmo/data/legacysurvey'
    return os.path.join(proj,drNumber)

def is_bool(obj):                                        
    return obj.dtype == bool

def is_numeric(obj):                                        
    try: 
        tmp=obj+5
    except TypeError:
        return False
    return True


def flux2mag(nmgy):
    return -2.5 * (np.log10(nmgy) - 9) 

class Bit(object):
    def set(self,value, bit):
        """change bit to 1, bit is 0-indexed"""
        return value | (1<<bit)

    def clear(self,value, bit):
        """change bit to 0, bit is 0-indexed"""
        return value & ~(1<<bit)


class TargetSelection(object):
    def __init__(self,prefix=''): 
        """Returns bool array of either eboss or desi ELGs
        
        Args:
            prefix: 'tractor_' for randoms table, '' for tractor table
        """
        self.prefix=prefix

    def run(self,tractor,name):
        assert(name in ['desi','eboss_ngc','eboss_sgc'])
        if name == 'desi':
            return self.desi_iselg(tractor)
        elif name == 'eboss_ngc':
            return self.eboss_iselg(tractor,'ngc')
        elif name == 'eboss_sgc':
            return self.eboss_iselg(tractor,'sgc')

    def desi_iselg(self,tractor):
        kw={}
        if self.prefix+'brick_primary' in tractor.get_columns():
            kw.update(primary=tractor.get(self.prefix+'brick_primary'))
        for band,iband in [('g',1),('r',2),('z',4)]:
            kw[band+'flux']= tractor.get(self.prefix+'flux_'+band) / \
                              tractor.get(self.prefix+'mw_transmission_'+band)
        return self._desi_iselg(**kw)

    def eboss_iselg(self,tractor,n_or_s):
        kw=dict(n_or_s=n_or_s)
        if self.prefix+'brick_primary' in tractor.get_columns():
            kw.update(primary=tractor.get(self.prefix+'brick_primary'))
        for key in ['ra','dec']:
            kw[key]= tractor.get(key) #self.prefix+key)
        for band in 'grz':
            kw['anymask_'+band]= tractor.get(self.prefix+'anymask_'+band)
        kw.update( self.get_grz_mag_dict(tractor) )
        return self._eboss_iselg(**kw)
 

    def _desi_iselg(self,gflux=None, rflux=None, zflux=None, 
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
    
    def _eboss_iselg(self,n_or_s=None,
                     primary=None,ra=None,dec=None,
                     gmag=None,rmag=None,zmag=None,
                     anymask_g=None,anymask_r=None,anymask_z=None):
        if primary is None:
            primary = np.ones(len(ra), bool)
        #print('dec.min',dec.min(),'dec.max',dec.max())
        #inRegion=dict(ngc= ((primary) &
        #                    (ra > 126.) &
        #                    (ra < 168.) &
        #                    (dec > 14.) &
        #                    (ra < 34.)),
        #              sgc_a= ((primary) &
        #                      (ra > 317.) &
        #                      (ra < 360.) &
        #                      (dec > -2.) &
        #                      (ra < 2.)),
        #              sgc_b= ((primary) &
        #                      (ra > 0.) &
        #                      (ra < 45.) &
        #                      (dec > -5.) &
        #                      (ra < 5.)))
        #inRegion.update(sgc= ((inRegion['sgc_a']) | 
        #                      (inRegion['sgc_b'])))
        # tycho2inblob == False
        # SDSS bright object mask & 0 < V < 11.5 mag Tycho2 stars mask
        # anymask[grz] == 0
        # custom mask for eboss23
        gr= gmag - rmag
        rz= rmag - zmag
        if n_or_s == 'ngc':
            colorCut= ((gmag > 21.825) &
                       (gmag < 22.9) &
                       (-0.068 * rz + 0.457 < gr) &
                       (gr < 0.112 * rz + 0.773) &
                       (0.637 * gr + 0.399 < rz) &
                       (rz < -0.555 * gr + 1.901))
        elif n_or_s == 'sgc':
            colorCut= ((gmag > 21.825) &
                       (gmag < 22.825) &
                       (-0.068 * rz + 0.457 < gr) &
                       (gr < 0.112 * rz + 0.773) &
                       (0.218 * gr + 0.571 < rz) &
                       (rz < -0.555 * gr + 1.901))
        anymask= ((anymask_g == 0) & 
                  (anymask_r == 0) & 
                  (anymask_z == 0))
 
        return (colorCut) & (anymask)

    def get_grz_mag_dict(self,tractor):
        d={}
        for band,iband in [('g',1),('r',2),('z',4)]:
            flux_ext= tractor.get(self.prefix+'flux_'+band) / \
                       tractor.get(self.prefix+'mw_transmission_'+band)
            d[band+'mag']= flux2mag(flux_ext)
        return d




class RandomsTable(object):
    """Creates the uniform,obiwan_a,obiwan_b randoms tables for a single brick

    Final table has same number rows as uniform table, and the obiwan
        rows are filled in wherever there is a matching unique_id
        between uniform and obiwan. A bitmask column 'obiwan_mask' which says
        whether the random was recovered by Tractor or not, and whether
        the random is near a previously existing real source in a DR 
        catalogue, like DR3 or DR5
    """
    def __init__(self, data_dir,dr3_or_dr5,db_randoms_table,
                 date='mm-dd-yyyy'):
        self.data_dir= data_dir
        self.dr3_or_dr5= dr3_or_dr5
        self.db_randoms_table= db_randoms_table
        self.date= date
        self.number_rsdirs= self.num_rsdirs_for_completion()

    def run(self,brick):
        tab= self.merge_randoms_tables(brick)
        if tab:
            self._run(tab,brick)
        else:
            print('brick %s skipped' % brick)

    def _run(self,tab,brick):
        self.add_flag_for_realsources(tab,brick)
        self.add_targets_mask(tab)
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
        if len(rsdirs) == self.number_rsdirs:
            return self._merge_randoms_tables(rsdirs)
        elif len(rsdirs) < self.number_rsdirs:
            print('brick %s not complete, %d/%d rsdirs exists' % \
                  (brick,len(rsdirs),self.number_rsdirs))
            return None
        else:
            raise ValueError('brick %s more rsdirs than should be possible %d/%d' % \
                             (brick,len(rsdirs),self.number_rsdirs))
    
    def _merge_randoms_tables(self,rsdirs):
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
            simcat= self.add_psql_to_uniform_table(simcat,self.db_randoms_table)
            # Recovered by Tractor
            tractor= fits_table(os.path.join(dr,'tractor-%s.fits' % brick))
            tractor.cut(tractor.brick_primary)
            cols= np.array(tractor.get_columns())
            del_cols= cols[(pd.Series(cols)
                              .str.startswith('apflux_'))]
            for col in del_cols:
                if col in ['apflux_resid_g','apflux_resid_r','apflux_resid_z',
                           'apflux_g','apflux_r','apflux_z']:
                    continue
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
                if is_bool(tractor.get(trac_key)):
                    add_vals[key]= np.zeros(len(simcat),bool)
                elif is_numeric(tractor.get(trac_key)):
                    shp= (len(simcat),) + tractor.get(trac_key).shape[1:]
                    add_vals[key]= np.zeros(shp) +np.nan
                else:
                    add_vals[key]= np.array([4*' ']*len(simcat))
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
        if any(db_dict['id'] - uniform.id != 0):
           sort_db= np.argsort(db_dict['id'])
           sort_uniform= np.argsort(uniform.id)
           for key in db_dict.keys():
               db_dict[key]= db_dict[key][sort_db]
           uniform= uniform[sort_uniform]
        assert(all(db_dict['id'] - uniform.id == 0))
        for key,val in db_dict.items():
            if key in ['id']:
                pass
            uniform.set('psql_%s' % key,val)
        return uniform

    def add_flag_for_realsources(self,tab,brick):
        """Flag sources also in DR3, DR5

        Args:
            tab: table returned by merged_randoms_table()
        """
        real= fits_table(os.path.join(datarelease_dir(self.dr3_or_dr5),
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

    def add_targets_mask(self,table):
        TS= TargetSelection(prefix='tractor_') 
        mask= np.zeros(len(table),dtype=np.int8)
        for survey_ts,bit in [('desi',0),
                              ('eboss_ngc',1),
                              ('eboss_sgc',2)]:
            keep= TS.run(table,survey_ts)
            if len(table[keep]) > 0:
                mask[keep]= Bit().set(mask[keep],bit)
        table.set('targets_mask',mask)

    def write_table(self,tab,fn):
        """Write the merged randoms table is doesn't already exist"""
        if not os.path.exists(fn): 
            tab.writeto(fn)
            print('Wrote %s' % fn)

    def num_rsdirs_for_completion(self):
        outdir_name= os.path.basename(self.data_dir)
        if outdir_name == 'eboss_elg':
            n=1
        elif outdir_name == 'elg_dr5_500per':
            n=7
        elif outdir_name == 'elg_dr5_1000per':
            n=4
        else:
            raise ValueError('%s not one of the above options' % outdir_name)
        return n

def fraction_recovered(randoms):
    return len(randoms[randoms.obiwan_mask == 1]) / float(len(randoms))

def bin_by_mag(randoms, func_to_apply, band=None,bin_minmax=(18.,26.),nbins=20):
    '''bins data and result of func_to_apply(randoms) into the bucket

    Args:
        randoms: randoms.fits table
        func_to_apply: operates on a randoms table, return val to store in bucket
        band: bin by mag in this band
    '''
    assert(band in 'grz')
    mag= flux2mag(randoms.get(band+'flux') / randoms.get('mw_transmission_'+band))

    bin_edges= np.linspace(bin_minmax[0],bin_minmax[1],num= nbins+1)
    vals={}
    vals['binc']= (bin_edges[1:]+bin_edges[:-1])/2.
    vals['val']=np.zeros(nbins)+np.nan
    for i,low,hi in zip(range(nbins), bin_edges[:-1],bin_edges[1:]):
        keep= ((mag > low) &
               (mag <= hi))
        if np.where(keep)[0].size > 0:
            vals['val'][i]= func_to_apply(randoms[keep])
        else:
            vals['val'][i]=np.nan 
    return vals



def depth_at_half_recovered(randoms,band):
    """bin by mag in given mag, return bin center where frac recovered first drops below 50%"""
    binc_and_vals= bin_by_mag(randoms, func_to_apply= fraction_recovered,
                              band=band,bin_minmax=(18.,26.),nbins=20)
    for i,val in enumerate(binc_and_vals['val']):
        if not np.isfinite(val):
            continue
        elif val <= 0.5:
            break
    #if i == len(binc_and_vals['val'])-1 and val > 0.5:
    if val > 0.5:
        raise ValueError('val=',val)
        return np.nan
    else:
        return binc_and_vals['binc'][i]

class SummaryTable(object):
    """Writes one table per brick, with brick averaged quantities
    
    derived table "randoms.fits" must exist. Joins the brick summary 
    quantities from a data release with a similar set from the 
    randoms.fits table. Each brick's table has one 
    row and all tables get merged to make the eatmap plots
    """
    def __init__(self, data_dir,dr3_or_dr5,date='mm-dd-yyyy'):
        self.data_dir= data_dir
        self.dr3_or_dr5= dr3_or_dr5
        self.date= date
        self.surveyBricks = fits_table(os.path.join(os.environ['LEGACY_SURVEY_DIR'],
                                                    'survey-bricks.fits.gz'))

    def run(self,brick):
        #summary_DR= self.brick_summary([brick])
        #summary_obi= self.brick_summary_obiwan(brick,prefix='tractor_')
        summary_obi= self.brick_summary_obiwan_brief(brick,prefix='tractor_')
        #self.add_obiwan_to_DR_table(summary_DR,summary_obi)
        # Write
        derived_dir= derived_field_dir(brick,self.data_dir,self.date)
        fn= os.path.join(derived_dir,'summary.fits')
        #self.write_table(summary_DR,fn)
        self.write_table(summary_obi,fn)

    def write_table(self,tab,fn):
        if not os.path.exists(fn): 
            tab.writeto(fn)
            print('Wrote %s' % fn)

    def add_obiwan_to_DR_table(self,summary_DR,summary_obi):
        """adds the summary_obi columsn to the summary_DR table
        
        Args:
            summary_DR: brick summary for the data release bricks
            summary_obiwan: brick summary for the obiwan bricks
        """
        prefix='obiwan_'
        for col in summary_obi.get_columns():
            summary_DR.set(prefix+col, summary_obi.get(col))
        del summary_obi

    def brick_summary_obiwan_brief(self,brick,prefix=''):
        """brick_summary_obiwan but only 3-4 quantities"""
        randoms_fn= os.path.join(derived_field_dir(brick,
                                        self.data_dir,self.date),
                                 'randoms.fits')
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
    
    def brick_summary_obiwan(self,brick,prefix=''):
        """brick summary for obiwan 
        
        Args:
            prefix: prefix for obiwan tractor columns, e.g. tractor_
        """
        randoms_fn= os.path.join(derived_field_dir(brick,
                                        self.data_dir,self.date),
                                 'randoms.fits')
        T = fits_table(randoms_fn)
        
        brickset = set()
        summary= defaultdict(list)
        nnhist = 6
    
        # Obiwan stuff
        was_recovered= T.obiwan_mask == 1
        summary['frac_recovered'].append( len(T[was_recovered])/ float(len(T)))
        for b in 'grz':
            summary['depth_at_half_recovered_'+b].append( depth_at_half_recovered(T,band=b))

        W = H = 3600
        # H=3600
        # xx,yy = np.meshgrid(np.arange(W), np.arange(H))
        unique = np.ones((H,W), bool)
        tlast = 0
        
        brickset.add(brick)
        for key in ['gn','rn','zn']:
            summary[key].append(0)

        for key in ['gnhist','rnhist','znhist']:
            summary[key].append([0 for i in range(nnhist)])

        ibrick = np.nonzero(self.surveyBricks.brickname == brick)[0][0]
        summary['ibricks'].append(ibrick)

        #T.cut(T.brick_primary)
        summary['nsrcs'].append(len(T))
        types = Counter([t.strip() for t in T.get(prefix+'type')])
        for typ in 'psf simp rex exp dev comp'.split(' '):
            summary['n'+typ].append(types[typ.upper()])
        print('N sources', summary['nsrcs'][-1])

        for b in 'grz':
            summary[b+'psfsize'].append(np.median(T.get(prefix+'psfsize_'+b)))
            summary[b+'psfdepth'].append(np.median(T.get(prefix+'psfdepth_'+b)))
            summary[b+'galdepth'].append(np.median(T.get(prefix+'galdepth_'+b)))
            summary[b+'trans'].append(np.median(T.get(prefix+'mw_transmission_'+b)))
            
        summary['ebv'].append(np.median(T.get(prefix+'ebv')))

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
         
        ibricks = np.array(summary['ibricks'])
        
        #print('Maximum number of sources:', max(nsrcs))
        
        T = fits_table()
        #T.brickname = np.array([brick])
        #T.ra  = self.surveyBricks.ra [ibricks]
        #T.dec = self.surveyBricks.dec[ibricks]
        T.set('frac_recovered', np.array(summary['frac_recovered']).astype(np.float32))
        for b in 'grz':
            T.set('depth_at_half_recovered_'+b, np.array(summary['depth_at_half_recovered_'+b]).astype(np.float32))
            T.set('nexp_'+b, np.array(summary[b+'n']).astype(np.int16))
            T.set('nexphist_'+b, np.array(summary[b+'nhist']).astype(np.int32))
            T.set('nobjs', np.array(summary['nsrcs']).astype(np.int16))
            T.set('psfsize_'+b, np.array(summary[b+'psfsize']).astype(np.float32))
            T.set('trans_'+b, np.array(summary[b+'trans']).astype(np.float32))
            T.set('ext_'+b, -2.5 * np.log10(T.get('trans_'+b)))
        
        for typ in 'psf simp rex exp dev comp'.split(' '):
            T.set('n'+typ, np.array(summary['n'+typ]).astype(np.int16))
        
        with np.errstate(divide='ignore'):
            for b in 'grz':
                T.set('psfdepth_'+b, (-2.5*(-9.+np.log10(5.*np.sqrt(1. / np.array(summary[b+'psfdepth']))))).astype(np.float32))
                T.set('galdepth_'+b, (-2.5*(-9.+np.log10(5.*np.sqrt(1. / np.array(summary[b+'galdepth']))))).astype(np.float32))
        
        for k in ['psfdepth_g', 'psfdepth_r', 'psfdepth_z', 
                  'galdepth_g', 'galdepth_r', 'galdepth_z']:
            v = T.get(k)
            v[np.logical_not(np.isfinite(v))] = 0.
        
        T.ebv = np.array(summary['ebv']).astype(np.float32)
        return T 




    def brick_summary(self,bricklist=[]):
        """
        Args:
            bricklist: Give a single brick as a list of length 1, e.g. [brick]
        """
        assert(len(bricklist) == 1)
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
       
        dirprefix= datarelease_dir(self.dr3_or_dr5)
        for ibrick,brick in enumerate(bricklist):
            #words = fn.split('/')
            #dirprefix = '/'.join(words[:-4])
            #print('Directory prefix:', dirprefix)
            #words = words[-4:]
            #brick = words[2]
            #print('Brick', brick)
            tfn = os.path.join(dirprefix, 'tractor', brick[:3], 'tractor-%s.fits'%brick)
            if self.dr3_or_dr5 == 'dr5': 
                columns=['brick_primary', 'type',
                         'psfsize_g', 'psfsize_r', 'psfsize_z',
                         'psfdepth_g', 'psfdepth_r', 'psfdepth_z',
                         'galdepth_g', 'galdepth_r', 'galdepth_z',
                         'ebv',
                         'mw_transmission_g', 'mw_transmission_r', 'mw_transmission_z',
                         'nobs_w1', 'nobs_w2', 'nobs_w3', 'nobs_w4',
                         'nobs_g', 'nobs_r', 'nobs_z',
                         'mw_transmission_w1', 'mw_transmission_w2', 'mw_transmission_w3', 'mw_transmission_w4']
            elif self.dr3_or_dr5 == 'dr3': 
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


            if self.dr3_or_dr5 == 'dr5':
                hasBands= [band for band in 'grz' if any(T.get('nobs_'+band) > 0)]
            elif self.dr3_or_dr5 == 'dr3':
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

            if self.dr3_or_dr5 == 'dr5':
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
                
            elif self.dr3_or_dr5 == 'dr3':
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

def main_mpi(bricks=[],doWhat=None,dr3_or_dr5=None,
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
        tabMaker= RandomsTable(data_dir,dr3_or_dr5,db_randoms_table,
                               date=date)
    elif doWhat == 'summary':
        tabMaker= SummaryTable(data_dir,dr3_or_dr5,date=date)
    
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
    parser.add_argument('--doWhat', type=str, choices=['randoms','summary'],required=True)
    parser.add_argument('--data_dir', type=str, required=True, 
                        help='path to obiwan/, tractor/ dirs') 
    parser.add_argument('--db_randoms_table', type=str, choices=['obiwan_eboss_elg',
                                    'obiwan_elg_dr5','obiwan_cosmos'],required=False)
    parser.add_argument('--nproc', type=int, default=1, help='set to > 1 to run mpi4py') 
    parser.add_argument('--bricks_fn', type=str, default=None,
                        help='specify a fn listing bricks to run, or a single default brick will be ran') 
    parser.add_argument('--dr3_or_dr5', type=str, choices=['dr5','dr3'], 
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

