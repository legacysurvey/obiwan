# See LICENSE.rst for BSD 3-clause license info
# -*- coding: utf-8 -*-
"""
====================
obiwan.priors
====================

Uses intrinsic (dust removed) AB fluxes/mags to select ELGs,LRGs,QSOs,STARs
Makes color-color plots that reproduce the FDR
Single band, mag distributions, plotted using 'as observed' AB fluxes/mags,
since these are what one adds into a CP image
"""

import matplotlib
if __name__ == '__main__':
    matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
from astropy.io import fits
#from astropy.table import vstack, Table
from astrometry.util.fits import fits_table, merge_tables
import os
import sys
from glob import glob
from scipy.optimize import newton
from sklearn.neighbors import KernelDensity
import pickle
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from six.moves import urllib
import tarfile
import pandas as pd
from pandas.plotting import scatter_matrix

from obiwan.common import inJupyter, save_png, fits2pandas
from obiwan.fetch import fetch_targz
from theValidator.catalogues import CatalogueFuncs,Matcher

DOWNLOAD_ROOT = "http://portal.nersc.gov/project/desi/users/kburleigh/obiwan/"
NERSC_ROOT = DOWNLOAD_ROOT.replace("http://portal.nersc.gov/project/",
                                   "/global/project/projectdirs/")\
                          .replace("/users/","/www/users/")
OBJTYPES= ['elg','lrg','star','qso']

plt.rcParams['axes.labelsize'] = 14
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12

class EmptyClass(object):
    pass

class Data(object):
    """Fetches and loads truth and catalogue data for star,galaxy populations
    
    Stored at: /global/project/projectdirs/desi/www/users/kburleigh/obiwan/priors/
    Originally at: dr2-dir, /project/projectdirs/desi/target/analysis/truth
    dr3-dir, /project/projectdirs/desi/users/burleigh/desi/target/analysis/truth
    
    Args:
      img_cuts: dict, set keys to True/False to turn imaging cuts On/Off
    """
    
    def __init__(self,**img_cuts):
	self.brick_primary= img_cuts.get('brick_primary',True)
	self.anymask= img_cuts.get('anymask',False)
	self.allmask= img_cuts.get('allmask',False)
	self.fracflux= img_cuts.get('fracflux',False)
        
    def fetch(self,outdir):
        name= 'priors.tar.gz'
        fetch_targz(os.path.join(DOWNLOAD_ROOT, name), outdir)
        
    def load_elg(self,DR, rlimit=23.4+1):
        """
        Args: 
            DR: which Data Release's data to use
            rlimit: AB magnitude limit of sample
                default: is 1 mag deeper than DESI requirement
        """
    	if DR == 2:
    	    zcat = fits_table(os.path.join(self.outdir,'deep2-field1-oii.fits.gz'))
    	    R_MAG= zcat.get('cfhtls_r')
    	    rmag_cut= R_MAG < rlimit 
    	elif DR == 3:
    	    zcat = fits_table(os.path.join(self.outdir,'deep2f234-dr3matched.fits'))
    	    decals = fits_table(os.path.join(self.outdir,'dr3-deep2f234matched.fits'))
            # Add mag data 
    	    CatalogueFuncs().set_mags_OldDataModel(decals)
    	    rflux= decals.get('decam_flux_nodust')[:,2]
    	    rmag_cut= rflux > 10**((22.5 - rlimit)/2.5)
    	    rmag_cut*= self.imaging_cut(decals)
    	zcat.cut(rmag_cut) 
    	decals.cut(rmag_cut) 
    	# color data
    	if DR == 2:
    	    G_MAG= zcat.get('cfhtls_g')
    	    R_MAG= zcat.get('cfhtls_r')
    	    Z_MAG= zcat.get('cfhtls_z')
    	elif DR == 3:
    	    G_MAG= decals.get('decam_mag_wdust')[:,1]
    	    R_MAG= decals.get('decam_mag_wdust')[:,2]
    	    Z_MAG= decals.get('decam_mag_wdust')[:,4]
    	    R_MAG_nodust= decals.get('decam_mag_nodust')[:,2]
    	    # Don't need w1, but might be useful
    	    W1_MAG = decals.get('wise_mag_wdust')[:,0]
    	# Repackage
    	tab= fits_table()
    	# Deep2
    	tab.set('ra', zcat.ra)
    	tab.set('dec', zcat.dec)
    	tab.set('zhelio', zcat.zhelio)
    	tab.set('oii_3727', zcat.oii_3727)
    	tab.set('oii_3727_err', zcat.oii_3727_err)
    	# Decals
    	tab.set('g_wdust', G_MAG)
    	tab.set('r_wdust', R_MAG)
    	tab.set('z_wdust', Z_MAG)
    	tab.set('w1_wdust', W1_MAG) 
    	tab.set('r_nodust', R_MAG_nodust)
    	# DR3 shape info
    	for key in ['type','shapeexp_r','shapedev_r']:
    	    tab.set(key, decals.get(key))
            # Add shape info from tractor cats
            dic= self.get_tractor_shapes(tab)
            tab.set('tractor_re', dic['re'])
            tab.set('tractor_n', dic['n'])
    	return tab        
    
     # def load_elg_acs():
     #     print('Matching acs to dr3,deep2')
     #     deep= fits_table(os.path.join(self.outdir,
     #                                   'deep2f1234_acsgcmatched.fits'))
     #     acs= fits_table(os.path.join(self.outdir,
     #                                  'acsgc_deep2f1234matched.fits'))
     #     # Remove ra,dec from acs, then merge
     #     for key in ['ra','dec']:
     #         acs.delete_column(key)
     #     return merge_tables([deep,acs], columns='fillzero')
        
    def load_lrg(self,DR,zlimit=20.46+1):
        """
        Args: 
            DR: which Data Release's data to use
            zlimit: AB magnitude limit of sample
                default: is 1 mag deeper than DESI requirement
        """
    	if DR == 2:
    	    decals= fits_table(os.path.join(self.outdir,'decals-dr2-cosmos-zphot.fits.gz') )
    	    spec=self.read_fits( os.path.join(self.outdir,'cosmos-zphot.fits.gz') )
    	elif self.DR == 3:
    	    decals= fits_table( os.path.join(self.outdir,'dr3-cosmoszphotmatched.fits') )
    	    spec= fits_table( os.path.join(self.outdir,'cosmos-zphot-dr3matched.fits') )
    	# DECaLS
    	CatalogueFuncs().set_mags_OldDataModel(decals)
    	Z_FLUX = decals.get('decam_flux_nodust')[:,4]
    	W1_FLUX = decals.get('wise_flux_nodust')[:,0]
    	# Cuts
    	# BUG!!
    	keep= np.all((Z_FLUX > 10**((22.5 - zlimit)/2.5),\
    		      W1_FLUX > 0.),axis=0)
    	keep *= self.imaging_cut(decals)
    	decals.cut(keep) 
    	spec.cut(keep)
    	# Repackage
    	tab= fits_table()
    	# Cosmos
    	tab.set('ra', spec.ra)
    	tab.set('dec', spec.dec)
    	tab.set('zp_gal', spec.zp_gal)
    	tab.set('type_zphotcomos', spec.type)
    	tab.set('mod_gal', spec.mod_gal)
    	# DR3
    	tab.set('g_wdust', decals.get('decam_mag_wdust')[:,1])
    	tab.set('r_wdust', decals.get('decam_mag_wdust')[:,2])
    	tab.set('z_wdust', decals.get('decam_mag_wdust')[:,4])
    	tab.set('w1_wdust', decals.get('wise_mag_wdust')[:,0])
    	tab.set('r_nodust', decals.get('decam_mag_nodust')[:,2])
        tab.set('z_nodust', decals.get('decam_mag_nodust')[:,4])
        # DR3 shape info
        for key in ['type','shapeexp_r','shapedev_r']:
            tab.set(key, decals.get(key))
        # Add tractor shapes
    	dic= self.get_tractor_shapes(tab)
    	tab.set('tractor_re', dic['re'])
    	tab.set('tractor_n', dic['n'])
        return tab

    # def load_lrg_acs(self,DR,zlimit):
    #     cosmos= self.load_lrg(DR,zlimit)
    #     acs= fits_table(os.path.join(self.outdir,
    #                                  'acsgc_cosmosmatched.fits'))
    #     # WARNING: not sure what this file is
    #     # cosmos_acsgcmatched.fits
    #     # Remove ra,dec from acs, then merge
    #     for key in ['ra','dec']:
    #         acs.delete_column(key)
    #     return merge_tables([cosmos,acs], columns='fillzero')
    
    # def load_lrg_vipers(self,DR,zlimit=20.46+1):
    #     """LRGs from VIPERS in CFHTLS W4 field (12 deg2)
    #     """
    #     if DR == 2:
    #         decals = fits_table(os.path.join(self.outdir,'decals-dr2-vipers-w4.fits.gz'))
    #         vip = fits_table(os.path.join(self.outdir,'vipers-w4.fits.gz'))
    #     elif DR == 3:
    #         decals= fits_table(os.path.join(self.outdir,'dr3-vipersw1w4matched.fits'))
    #         vip = fits_table(os.path.join(self.outdir,'vipersw1w4-dr3matched.fits'))
    #     CatalogueFuncs().set_mags_OldDataModel(decals)
    #     Z_FLUX = decals.get('decam_flux_nodust')[:,4]
    #     W1_FLUX = decals.get('wise_flux_nodust')[:,0]
    #     index={}
    #     index['decals']= np.all((Z_FLUX > 10**((22.5 - zlimit)/2.5),\
    #     			 W1_FLUX > 0.),axis=0)
    #     index['decals']*= self.imaging_cut(decals)
    #     # VIPERS
    #     # https://arxiv.org/abs/1310.1008
    #     # https://arxiv.org/abs/1303.2623
    #     flag= vip.get('zflg').astype(int)
    #     index['good_z']= np.all((flag >= 2,\
    #     			 flag <= 9,\
    #     			 vip.get('zspec') < 9.9),axis=0) 
    #     # get Mags
    #     rz,rW1={},{}
    #     cut= np.all((index['decals'],\
    #     	     index['good_z']),axis=0)
    #     rz= decals.get('decam_mag_nodust')[:,2][cut] - decals.get('decam_mag_nodust')[:,4][cut]
    #     rW1= decals.get('decam_mag_nodust')[:,2][cut] - decals.get('wise_mag_nodust')[:,0][cut]
    #     return rz,rW1
    
    def load_star(self,DR):
    	if DR == 2:
    	    stars= fits_table(os.path.join(self.outdir,'Stars_str82_355_4.DECaLS.dr2.fits'))
    	    CatalogueFuncs().set_mags_OldDataModel(stars)
    	    stars.cut( self.std_star_cut(stars) )
    	elif self.DR == 3:
    	    raise ValueError()
            return stars

    # def load_star_sweep(self):
    #     """Model the g-r, r-z color-color sequence for stars"""
    #     # Build a sample of stars with good photometry from a single sweep.
    #     rbright = 18
    #     rfaint = 19.5
    #     swp_dir='/global/project/projectdirs/cosmo/data/legacysurvey/dr2/sweep/2.0'
    #     sweep = self.read_fits(os.path.join(swp_dir,'sweep-340p000-350p005.fits'))
    #     keep = np.where((sweep.get('type') == 'PSF ')*
    #     		(np.sum((sweep.get('decam_flux')[:, [1,2,4]] > 0)*1, axis=1)==3)*
    #     		(np.sum((sweep.get('DECAM_ANYMASK')[:, [1,2,4]] > 0)*1, axis=1)==0)*
    #     	        (np.sum((sweep.get('DECAM_FRACFLUX')[:, [1,2,4]] < 0.05)*1, axis=1)==3)*
    #     		(sweep.get('decam_flux')[:,2]<(10**(0.4*(22.5-rbright))))*
    #     		(sweep.get('decam_flux')[:,2]>(10**(0.4*(22.5-rfaint)))))[0]
    #     stars = sweep[keep]
    #     print('dr2stars sample: {}'.format(len(stars)))
    #     gg = 22.5-2.5*np.log10(stars.get('decam_flux')[:, 1])
    #     rr = 22.5-2.5*np.log10(stars.get('decam_flux')[:, 2])
    #     zz = 22.5-2.5*np.log10(stars.get('decam_flux')[:, 4])
    #     gr = gg - rr
    #     rz = rr - zz
    #	return np.array(rz),np.array(gr)

    
    def load_qso(self,DR,rlimit=22.7+1):
        """
        Args: 
            DR: which Data Release's data to use
            rlimit: AB magnitude limit of sample
                default: is 1 mag deeper than DESI requirement
        """
        if DR == 2:        
    	    qsos= fits_table( os.path.join(self.outdir,'AllQSO.DECaLS.dr2.fits') )
    	    # Add AB mags
    	    CatalogueFuncs().set_mags_OldDataModel(qsos)
    	    qsos.cut( self.imaging_cut(qsos) )
    	    # r < 22.7, grz > 17
    	    GFLUX = qsos.get('decam_flux_nodust')[:,1] 
    	    RFLUX = qsos.get('decam_flux_nodust')[:,2]
    	    ZFLUX = qsos.get('decam_flux_nodust')[:,4] 
    	    GRZFLUX = (GFLUX + 0.8* RFLUX + 0.5* ZFLUX ) / 2.3
    	    cut= np.all((RFLUX > 10**((22.5 - rlimit)/2.5),\
    			 GRZFLUX < 10**((22.5-17.0)/2.5)),axis=0)
    	    qsos.cut(cut)
    	elif DR == 3:
    	    raise ValueError('Not done yet')
            return qsos
    
    def imaging_cut(self,data):
    	"""data: tractor catalogue"""
    	cut=np.ones(len(data)).astype(bool)
    	# Brick Primary
    	if data.get('brick_primary').dtype == 'bool':
    	    cut*= data.get('brick_primary') == True
    	elif data.get('brick_primary').dtype == 'S1':
    	    cut*= data.get('brick_primary') == 'T'
    	else: 
    	    raise ValueError('brick_primary has type=',data.get('brick_primary').dtype)
    	#if self.anymask:
            #cut*= np.all((data.get('decam_anymask')[:, [1,2,4]] == 0),axis=1)
            # ALL Mask
            cut*= np.all((data.get('decam_allmask')[:, [1,2,4]] == 0),axis=1)
    	# FracFlux
    	cut*= np.all((data.get('decam_fracflux')[:, [1,2,4]] < 0.05),axis=1)
    	return cut

    def std_star_cut(self,data):
    	"""See: https://desi.lbl.gov/trac/wiki/TargetSelectionWG/TargetSelection#SpectrophotometricStandardStarsFSTD
    	data is a fits_table object with Tractor Catalogue columns
    	"""
    	RFLUX_obs = data.get('decam_flux')[:,2]
    	GFLUX = data.get('decam_flux_nodust')[:,1]
    	RFLUX = data.get('decam_flux_nodust')[:,2]
    	ZFLUX = data.get('decam_flux_nodust')[:,4]
    	GRZSN = data.get('decam_flux')[:,[1,2,4]] * np.sqrt(data.get('decam_flux_ivar')[:,[1,2,4]])
    	GRCOLOR = 2.5 * np.log10(RFLUX / GFLUX)
    	RZCOLOR = 2.5 * np.log10(ZFLUX / RFLUX)
    	cut= np.all((data.get('brick_primary') == True,\
    		     data.get('type') == 'PSF',\
    		     np.all((data.get('decam_allmask')[:, [1,2,4]] == 0),axis=1),\
    		     np.all((data.get('decam_fracflux')[:, [1,2,4]] < 0.04),axis=1),\
    		     np.all((GRZSN > 10),axis=1),\
    		     RFLUX_obs < 10**((22.5-16.0)/2.5)),axis=0)
    	#np.power(GRCOLOR - 0.32,2) + np.power(RZCOLOR - 0.13,2) < 0.06**2,\
    	#RFLUX_obs > 10**((22.5-19.0)/2.5)),axis=0)
    	return cut 
    
    def get_tractor_shapes(self,cat):
        """Get a dictionary of Tractor Catalogue shape measurements

        Args: 
            cat: tractor catalogues

        Returns: 
            d: dict of 
                {'re','n','ba','pa'}
        """
        d= {}
        for key in ['re','n','ba','pa']:
            d[key]= np.zeros(len(cat))-1
        # Re,N
        # SIMP
        #keep= (cat.type == 'SIMP') * (cat.shapeexp_r > 0.)
        #d['re'][keep]= cat.shapeexp_r[keep]
        #d['n'][keep]= 1.
        # EXP
        keep= (cat.type == 'EXP') * (cat.shapeexp_r > 0.)
        d['re'][keep]= cat.shapeexp_r[keep]
        d['n'][keep]= 1.
        # DEV
        keep= (cat.type == 'DEV') * (cat.shapedev_r > 0.)
        d['re'][keep]= cat.shapedev_r[keep]
        d['n'][keep]= 4.
        # BA, PA --> assume radmon dist, so don't return these
        #d['ba']= np.random.uniform(0.2,1.,size=len(cat))
        #d['pa']= np.random.uniform(0.,180.,size=len(cat))
        return d

def pyVersion():
    return sys.version_info.major
        
def save_pickle(obj,outdir,name):
    """Saves data in "obj" to pickle file"""
    # py2 and 3 pickle files not backward compatible
    path= os.path.join(outdir,name + "_py%d.pkl" % pyVersion())
    fout=open(path,'w')
    pickle.dump(obj,fout)
    fout.close()
    print("Saved pickle", path)

def load_pickle(outdir,name):
    """Loads data in pickle file"""
    # py2 and 3 pickle files not backward compatible
    path= os.path.join(outdir,name + "_py%d.pkl" % pyVersion())
    with open(path,'r') as fout:
        obj= pickle.load(fout)
    print("Loaded pickle", path)
    return obj    

 

class KDE_Model(object):
    """Fits a Kernel Density Estimate (KDE) to a given source population
    
    Only one source or data table per KDE_Model object 

    Args:
        src: elg,lrg, etc
        data: tractor catalogue for the population 
            returned by Data().load_elg(), for elgs 
        outdir: output dir
    
    Attributes:
        data: copy of data
        outdir:
        kde: fit KDE
    """
    
    def __init__(self,src,data,outdir):
        assert(src in OBJTYPES)
        self.src= src
        # DF instead of fits_table
        self.df= fits2pandas(data, attrs=self.keys_for_KDE(src))
        self.outdir= os.path.join(outdir,'kdes')
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

    def keys_for_KDE(self,src):
        assert src in OBJTYPES
        if src == "elg":
            return ['zhelio','oii_3727','oii_3727_err',
                    'r_wdust','z_wdust','g_wdust','tractor_re']
        elif src == 'lrg':
            return []
        elif src == 'star':
            return []
        elif src == 'qso':
            return []
        
    def get_kde(self):
        name= '%s-kde' % self.src
        try:
            self.df_wcut,self.df_for_kde,kde= load_pickle(self.outdir,name)
        except IOError:
            # Doesn't exist yet
            if self.src == 'elg':
                kde= self.fit_elg()
            elif self.src == 'lrg':
                kde= self.fit_lrg()
            elif self.src == 'star':
                kde= self.fit_star()
            elif self.src == 'qso':
                kde= self.fit_qso()
            save_pickle((self.df_wcut,self.df_for_kde,kde),self.outdir,name)
        return kde
    
    def fit_kde(self,X,
                bandwidth=0.05, kernel='tophat'):
        """Fits a KDE to the data
        
        Following 
            http://scikit-learn.org/stable/auto_examples/neighbors/plot_kde_1d.html
        
        Args: 
            X: numpy array of shape [m instances, N features]
            bandwidth: characteristic width of the kernel
            kernel: gaussian, tophat, etc
          
        Returns: 
            KernelDensity() object fit to data_list
        """
        assert(X.shape[0] > X.shape[1])
        return KernelDensity(kernel=kernel,
                             bandwidth=bandwidth).fit(X)
    
    def fit_elg(self):
        # CUTS
        # Sample that represents population
        print('fitting elg!!'.upper())
        keep= ((self.df.zhelio >= 0.8) &
	           (self.df.zhelio <= 1.4) &
	           (self.df.oii_3727 >= 0.) &
	           (self.df.oii_3727_err > 0.))
    	# Finite values
    	for col in self.df.keys():
    	    if col == 'type':
    		  continue
    	    keep *= (np.isfinite(self.df[col]))
        # Real size
    	keep*= (self.df.tractor_re > 0.)
    	self.df_wcut= self.df[keep]
        # DF to fit on
        d= dict(r_wdust=self.df_wcut.r_wdust,
                rz_wdust= self.df_wcut.r_wdust - self.df_wcut.z_wdust,
		gr_wdust= self.df_wcut.g_wdust - self.df_wcut.r_wdust,
                zhelio= self.df_wcut.zhelio,
		tractor_re= self.df_wcut.tractor_re)
        self.df_for_kde= pd.DataFrame(d)
        return self.fit_kde(self.df_for_kde.values,
			                bandwidth=0.05,kernel='tophat')

    def fit_lrg(self):
    	# Cut to sample that repr. population
    	self.data.cut(((self.data.type_zphotcomos == 0) &
    		           (self.data.mod_gal <= 8))) 
    	print('dr3_cosmos after obiwan cuts: %d' % len(self.data))
        # Only finite vals
    	keep= np.ones(len(self.data),bool) 
    	for name in self.data.get_columns():
    	    if name == 'type':
    		  continue
    	    keep *= (np.isfinite( self.data.get(name) ))
    	self.data.cut(keep)
    	print('dr3_cosmos after finite cuts: %d' % len(self.data))
    	# physical
    	tab.cut(tab.tractor_re > 0.)
    	# Reshift needs be > 0
    	bandwidth=0.05
    	self.data.cut(self.data.zp_gal - bandwidth >= 0.)
    	print('dr3_cosmos after redshift > %f: %d' % (bandwidth,len(self.data)))
    	return self.fit_kde([tab.z_wdust, tab.r_wdust - tab.z_wdust, 
    			     tab.r_wdust - tab.w1_wdust, tab.zp_gal, 
    			     tab.g_wdust, tab.tractor_re],
    			    bandwidth=bandwidth,kernel='tophat')

    def fit_star(self):
        return self.fit_kde([self.data.get('decam_mag_wdust')[:,2],
                             stars.get('decam_mag_nodust')[:,2]-stars.get('decam_mag_nodust')[:,4],
                             stars.get('decam_mag_nodust')[:,1]-stars.get('decam_mag_nodust')[:,2]],
			    bandwidth=0.05,kernel='gaussian')

    def fit_qso(self):
        # Sample
        hiz=2.1
        self.data.cut( (self.data.z <= hiz)*\
                       (self.data.z >= 0.))
        # Finite vals
        keep= np.ones(len(self.data),bool) 
    	for name in self.data.get_columns():
    	    if name == 'type':
    		  continue
    	    keep *= (np.isfinite( self.data.get(name) ))
    	self.data.cut(keep)
    	return self.fit_kde([self.data.get('decam_mag_wdust')[:,2],
                                 self.data.get('decam_mag_wdust')[:,2]-qsos.get('decam_mag_wdust')[:,4],
                                 self.data.get('decam_mag_wdust')[:,1]-qsos.get('decam_mag_wdust')[:,2],
                                 self.data.z],
                                bandwidth=0.05,kernel='gaussian') 

class KernelOfTruth(object):
    # 1D histograms of each dim of KDE
    def plot_indiv_1d(self,lims=None, ndraws=1000,prefix=''):
        samp= self.kde.sample(n_samples=ndraws)
        for i,name in enumerate(self.labels):
            fig,ax= plt.subplots()
            # Data
            h,edges= np.histogram(self.X[:,i],bins=40,normed=True)
            binc= (edges[1:]+edges[:-1])/2.
            ax.step(binc,h,where='mid',lw=1,c='k',label='Data')
            # KDE distribution
            h,edges= np.histogram(samp[:,i],bins=40,normed=True)
            binc= (edges[1:]+edges[:-1])/2.
            ax.step(binc,h,where='mid',lw=1,c='b',label='KDE')
            xlab=ax.set_xlabel(name,fontsize='x-large')
            ylab=ax.set_ylabel('PDF')
            if lims:
                ax.set_xlim(lims[i])
            ax.legend(loc='upper right')
            savenm= 'kde_1d_%s_%s.png' % (prefix,name)
            plt.savefig(savenm,bbox_extra_artists=[xlab,ylab], bbox_inches='tight',dpi=150)
            plt.close()
            print('Wrote %s' % savenm)

    # 2D scatterplot of selected dims of KDE
    def plot_FDR_using_kde(self,obj='LRG',ndraws=1000,prefix='',outdir=None,nb=False):
        assert(obj in ['LRG','ELG'])
        for use_data in [True,False]:
            fig,ax = plt.subplots()
            # Add box
            ts= TSBox(src=obj)
            xrange,yrange= xyrange['x_%s' % obj.lower()],xyrange['y_%s' % obj.lower()]
            ts.add_ts_box(ax, xlim=xrange,ylim=yrange)
            # xyrange dashed box
            #for i in range(2):
            #    ax[i].plot([xrange[i],xrange[i]],yrange,'k--')
            #    ax[i].plot(xrange,[yrange[i],yrange[i]],'k--')
            # KDE sample
            if obj == 'LRG':
                xname= 'rz'
                yname= 'rw1'
            elif obj == 'ELG':
                xname= 'rz'
                yname= 'gr'
            ix= np.where(self.labels == xname)[0][0]
            iy= np.where(self.labels == yname)[0][0]
            if use_data:
                ax.scatter(self.X[:,ix],self.X[:,iy],
                           c='k',marker='o',s=10.,rasterized=True,label='Data')
                savenm= 'kde_in_FDR_Data_%s_%s.png' % (prefix,obj)
            else:
                samp= self.kde.sample(n_samples=ndraws)
                ax.scatter(samp[:,ix],samp[:,iy],
                           c='b',marker='o',s=10.,rasterized=True,label='KDE')
                savenm= 'kde_in_FDR_KDE_%s_%s.png' % (prefix,obj)
            # finish
            #ax.set_xlim(xrange[0]-1,xrange[1]+1)
            #ax.set_ylim(yrange[0]-1,yrange[1]+1)
            ax.set_xlim(xrange)
            ax.set_ylim(yrange)
            ax.legend(loc='upper right')
            if obj == 'LRG':
                xlab=ax.set_xlabel('r-z')
                ylab=ax.set_ylabel('r-W1')
            elif obj == 'ELG':
                xlab=ax.set_xlabel('r-z')
                ylab=ax.set_ylabel('g-r')
            ax.set_aspect(1)
            if outdir:
                savenm= os.path.join(outdir,savenm)
            plt.savefig(savenm,bbox_extra_artists=[xlab,ylab], bbox_inches='tight',dpi=150)
            print('Wrote %s' % savenm)
            if not nb:
                plt.close()

    def plot_1band_and_color(self, ndraws=1000,xylims=None,prefix=''):
        """xylims -- dict of x1,y1,x2,y2,... where x1 is tuple of low,hi for first plot xaxis"""
        fig,ax= plt.subplots(2,2,figsize=(15,10))
        plt.subplots_adjust(wspace=0.2,hspace=0.2)
        # Data
        xyz= self.kde.sample(n_samples=ndraws)
        ax[0,0].hist(self.X[:,0],normed=True)
        ax[0,1].scatter(self.X[:,1],self.X[:,2],\
                      c='b',edgecolors='none',marker='o',s=10.,rasterized=True,alpha=0.2)
        # KDE distribution
        xyz= self.kde.sample(n_samples=ndraws)
        ax[1,0].hist(xyz[:,0],normed=True)
        ax[1,1].scatter(xyz[:,1],xyz[:,2],\
                      c='b',edgecolors='none',marker='o',s=10.,rasterized=True,alpha=0.2)
        for cnt in range(2):
            if xylims is not None:
                ax[cnt,0].set_xlim(xylims['x1'])
                ax[cnt,0].set_ylim(xylims['y1'])
                ax[cnt,1].set_xlim(xylims['x2'])
                ax[cnt,1].set_ylim(xylims['y2'])
            xlab=ax[cnt,0].set_xlabel(self.labels[0],fontsize='x-large')
            xlab=ax[cnt,1].set_xlabel(self.labels[1],fontsize='x-large')
            ylab=ax[cnt,1].set_ylabel(self.labels[2],fontsize='x-large')
        plt.savefig('%skde.png' % prefix,bbox_extra_artists=[xlab], bbox_inches='tight',dpi=150)
        plt.close()
        if prefix == 'lrg_':
            # plot g band distribution even though no Targeting cuts on g
            fig,ax= plt.subplots(2,1,figsize=(8,10))
            plt.subplots_adjust(hspace=0.2)
            ax[0].hist(self.X[:,3],normed=True)
            ax[1].hist(xyz[:,3],normed=True) 
            for cnt in range(2):
                if xylims is not None:
                    ax[cnt].set_xlim(xylims['x3'])
                    ax[cnt].set_ylim(xylims['y3'])
            xlab=ax[0].set_xlabel(self.labels[3],fontsize='x-large')
            plt.savefig('%sg_kde.png' % prefix,bbox_extra_artists=[xlab], bbox_inches='tight',dpi=150)
            plt.close()

    def plot_1band_color_and_redshift(self, ndraws=1000,xylims=None,prefix=''):
        """xylims -- dict of x1,y1,x2,y2,... where x1 is tuple of low,hi for first plot xaxis"""
        fig,ax= plt.subplots(2,3,figsize=(20,10))
        plt.subplots_adjust(wspace=0.2,hspace=0.2)
        # Colormap the color-color plot by redshift
        cmap = mpl.colors.ListedColormap(['m','r', 'y', 'g','b', 'c'])
        bounds= np.linspace(xylims['x3'][0],xylims['x3'][1],num=6)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        # Data
        ax[0,0].hist(self.X[:,0],normed=True)
        # color bar with color plot
        axobj= ax[0,1].scatter(self.X[:,1],self.X[:,2],c=self.X[:,3],\
                                 marker='o',s=10.,rasterized=True,lw=0,\
                                 cmap=cmap,norm=norm,\
                                 vmin=bounds.min(),vmax=bounds.max())
        divider3 = make_axes_locatable(ax[0,1])
        cax3 = divider3.append_axes("right", size="5%", pad=0.1)
        cbar3 = plt.colorbar(axobj, cax=cax3,\
                             cmap=cmap, norm=norm, boundaries=bounds, ticks=bounds)
        cbar3.set_label('redshift')
        #
        ax[0,2].hist(self.X[:,3],normed=True)
        #Xtest= np.linspace(-0.2,1.6,num=100)
        #log_dens = self.kde.score_samples(self.X)
        #ax[0,2].plot(self.X[:,3], np.exp(log_dens), 'y-')
        # KDE distribution
        samp= self.kde.sample(n_samples=ndraws)
        ax[1,0].hist(samp[:,0],normed=True)
        ax[1,1].scatter(samp[:,1],samp[:,2],\
                      c='b',edgecolors='none',marker='o',s=10.,rasterized=True,alpha=0.2)
        ax[1,2].hist(samp[:,3],normed=True)
        for cnt in range(2):
            if xylims is not None:
                ax[cnt,0].set_xlim(xylims['x1'])
                ax[cnt,0].set_ylim(xylims['y1'])
                ax[cnt,1].set_xlim(xylims['x2'])
                ax[cnt,1].set_ylim(xylims['y2'])
                ax[cnt,2].set_xlim(xylims['x3'])
                ax[cnt,2].set_ylim(xylims['y3'])
            xlab=ax[cnt,0].set_xlabel(self.labels[0],fontsize='x-large')
            xlab=ax[cnt,1].set_xlabel(self.labels[1],fontsize='x-large')
            ylab=ax[cnt,1].set_ylabel(self.labels[2],fontsize='x-large')
            xlab=ax[cnt,2].set_xlabel(self.labels[3],fontsize='x-large')
        sname='%skde.png' % prefix
        plt.savefig(sname,bbox_extra_artists=[xlab], bbox_inches='tight',dpi=150)
        plt.close()
        print('Wrote %s' % sname)
        if prefix == 'lrg_':
            # plot g band distribution even though no Targeting cuts on g
            fig,ax= plt.subplots(2,1,figsize=(8,10))
            plt.subplots_adjust(hspace=0.2)
            ax[0].hist(self.X[:,4],normed=True)
            ax[1].hist(samp[:,4],normed=True) 
            for cnt in range(2):
                if xylims is not None:
                    ax[cnt].set_xlim(xylims['x4'])
                    ax[cnt].set_ylim(xylims['y4'])
            xlab=ax[0].set_xlabel(self.labels[4],fontsize='x-large')
            sname='%sg_kde.png' % prefix
            plt.savefig(sname,bbox_extra_artists=[xlab], bbox_inches='tight',dpi=150)
            plt.close()
            print('Wrote %s' % sname)

    def plot_galaxy_shapes(self, ndraws=1000,xylims=None,name='kde.png'):
        """xylims -- dict of x1,y1,x2,y2,... where x1 is tuple of low,hi for first plot xaxis"""
        fig,ax= plt.subplots(2,4,figsize=(20,10))
        plt.subplots_adjust(wspace=0.2,hspace=0.2)
        samp= self.kde.sample(n_samples=ndraws)
        # ba,pa can be slightly greater 1.,180
        samp[:,2][ samp[:,2] > 1. ]= 1.
        samp[:,3][ samp[:,3] > 180. ]= 180.
        # Physical values
        assert(np.all(samp[:,0] > 0))
        assert(np.all((samp[:,1] > 0)*\
                      (samp[:,1] < 10)))
        assert(np.all((samp[:,2] > 0)*\
                      (samp[:,2] <= 1.)))
        assert(np.all((samp[:,3] >= 0)*\
                      (samp[:,3] <= 180)))
        # plot
        for cnt in range(4):
            if cnt == 0:
                bins=np.linspace(0,80,num=20)
                # Data
                ax[0,cnt].hist(self.X[:,cnt],bins=bins,normed=True)
                # KDE distribution
                ax[1,cnt].hist(samp[:,cnt],bins=bins,normed=True)
            else:
                ax[0,cnt].hist(self.X[:,cnt],normed=True)
                ax[1,cnt].hist(samp[:,cnt],normed=True)
        # lims
        for row in range(2):
            for col in range(4):
                if xylims is not None:
                    ax[row,col].set_xlim(xylims['x%s' % str(col+1)])
                    #ax[cnt,1].set_xlim(xylims['x2'])
                    #ax[cnt,2].set_xlim(xylims['x3'])
                    xlab=ax[row,col].set_xlabel(self.labels[col],fontsize='x-large')
                    #xlab=ax[cnt,1].set_xlabel(self.labels[1],fontsize='x-large')
        plt.savefig(name,bbox_extra_artists=[xlab], bbox_inches='tight',dpi=150)
        plt.close()
        print('Wrote %s' % name)



def get_rgb_cols():
    return [(255,0,255),(102,255,255),(0,153,153),\
            (255,0,0),(0,255,0),(0,0,255),\
            (0,0,0)]


# Globals
xyrange=dict(x_star=[-0.5,2.2],\
             y_star=[-0.3,2.],\
             x_elg=[-0.5,2.2],\
             y_elg=[-0.3,2.],\
             x_lrg= [0, 3.],\
             y_lrg= [-2, 6],\
             x1_qso= [-0.5,3.],\
             y1_qso= [-0.5,2.5],\
             x2_qso= [-0.5,4.5],\
             y2_qso= [-2.5,3.5])

def rm_last_ticklabel(ax):
    """for multiplot"""
    labels=ax.get_xticks().tolist()
    labels=np.array(labels).astype(float) #prevent from making float
    labels=list(labels)
    labels[-1]=''
    ax.set_xticklabels(labels)



        
    
class Plot(object):
    """Makes relavent plots for a given src using data,kde      
    
    Args:
        data: tractor catalogue for source population
        kde: Fit KDE returned by FitKDE()
    
    """
    def __init__(self,outdir):
        self.outdir= os.path.join(outdir,'plots')
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
            
    def elg(self,data,df_wcut,df_for_kde,kde):
        src= 'elg'
        # Pandas built in
        df_wcut.hist(bins=20,figsize=(12,8))
        save_png(self.outdir,'%s_hist_df_wcut' % src)
        scatter_matrix(df_wcut, figsize=(20, 15))
        save_png(self.outdir,'%s_scat_df_wcut' % src)
        ####
        df_for_kde.hist(bins=20,figsize=(12,8))
        save_png(self.outdir,'%s_hist_df_for_kde' % src)
        scatter_matrix(df_for_kde, figsize=(20, 15))
        save_png(self.outdir,'%s_scat_df_for_kde' % src)
        ###
        samp= kde.sample(n_samples=1000)
        df_samp= pd.DataFrame(samp,columns=(df_for_kde.keys()))
        df_samp.hist(bins=20,figsize=(20,15))
        save_png(self.outdir,'%s_hist_df_samp' % src)
        scatter_matrix(df_samp, figsize=(20, 15))
        save_png(self.outdir,'%s_scat_df_samp' % src)

    def lrg(self,data,kde):
        src='lrg'
        pass
    
    def star(self,data,kde):
        src='star'
        pass
    
    def qso(self,data,kde):
        src='qso'
        pass



class TSBox(object):
    """Add FDR Target Selection box to any  ELG, LRG, plot
    
	Args:
		src: ELG, LRG, QSO, STAR

	Attributes:
		src: ELG, LRG, QSO, STAR
		
	Note: 
		add_ts_box(): primary function to use
	"""
	
    def __init__(self,src='ELG'):
        self.src=src

    def add_ts_box(self, ax, xlim=None,ylim=None):
        """draw color selection box"""
        assert(xlim is not None and ylim is not None)
        if self.src == 'ELG':
            #g-r vs. r-z
            xint= newton(self.ts_root,np.array([1.]),args=('y1-y2',))
            
            x=np.linspace(xlim[0],xlim[1],num=1000)
            y=np.linspace(ylim[0],ylim[1],num=1000)
            x1,y1= x,self.ts_box(x,'y1')
            x2,y2= x,self.ts_box(x,'y2')
            x3,y3= np.array([0.3]*len(x)),y
            x4,y4= np.array([0.6]*len(x)),y
            b= np.all((x >= 0.3,x <= xint),axis=0)
            x1,y1= x1[b],y1[b]
            b= np.all((x >= xint,x <= 1.6),axis=0)
            x2,y2= x2[b],y2[b]
            b= y3 <= np.min(y1)
            x3,y3= x3[b],y3[b]
            b= y4 <= np.min(y2)
            x4,y4= x4[b],y4[b]
            ax.plot(x1,y1,'k--',lw=2)
            ax.plot(x2,y2,'k--',lw=2)
            ax.plot(x3,y3,'k--',lw=2)
            ax.plot(x4,y4,'k--',lw=2) 
        elif self.src == 'LRG':
            #r-w1 vs. r-z
            x=np.linspace(xlim[0],xlim[1],num=1000)
            y=np.linspace(ylim[0],ylim[1],num=1000)
            x1,y1= x,self.ts_box(x,'y1')
            x2,y2= np.array([1.5]*len(x)),y
            b= x >= 1.5
            x1,y1= x1[b],y1[b]
            b= y2 >= np.min(y1)
            x2,y2= x2[b],y2[b]
            ax.plot(x1,y1,'k--',lw=2)
            ax.plot(x2,y2,'k--',lw=2)
        else: raise ValueError('src=%s not supported' % src)

    def ts_box(self, x,name):
        if self.src == 'ELG':
            if name == 'y1': return 1.15*x-0.15
            elif name == 'y2': return -1.2*x+1.6
            else: raise ValueError
        elif self.src == 'LRG':
            if name == 'y1': return 1.8*x-1.
            else: raise ValueError
        else: raise ValueError('src=%s not supported' % self.src)

    def ts_root(self,x,name):
        if self.src == 'ELG':
            if name == 'y1-y2': return self.ts_box(x,'y1')-self.ts_box(x,'y2')
            else: raise ValueError
        else: raise ValueError('non ELG not supported')

        
class FDR_elg(object):    
    def get_FDR_cuts(self,tab):
    	oiicut1 = 8E-17 # [erg/s/cm2]
    	zmin = 0.6
    	keep= {}
    	keep['lowz'] = tab.zhelio < zmin
    	keep['medz_lowO2'] = np.all((tab.zhelio > zmin,\
    				     tab.oii_3727_err != -2.0,\
    				     tab.oii_3727 < oiicut1), axis=0)
    	keep['medz_hiO2'] = np.all((tab.zhelio > zmin,\
    				    tab.zhelio < 1.0,\
    				    tab.oii_3727_err != -2.0,\
    				    tab.oii_3727 > oiicut1), axis=0)
    	keep['hiz_hiO2'] = np.all((tab.zhelio > 1.0,\
    				   tab.oii_3727_err !=-2.0,\
    				   tab.oii_3727 > oiicut1),axis=0)
    	return keep 		 
    
        
    def plot_FDR(self):
    	tab= self.get_dr3_deep2()
    	# Cuts 'lowz','medz_lowO2' ...
    	keep= self.get_FDR_cuts(tab)
    	# Plot
    	fig, ax = plt.subplots()
    	# Add box
    	ts= TSBox(src='ELG')
    	xrange,yrange= xyrange['x_elg'],xyrange['y_elg']
        ts.add_ts_box(ax, xlim=xrange,ylim=yrange)
    	# Add points
    	for cut_name,lab,color,marker in zip(['lowz','medz_lowO2','medz_hiO2','hiz_hiO2'],
    					    [r'$z<0.6$',r'$z>0.6, [OII]<8\times10^{-17}$',
    					     r'$z>0.6, [OII]>8\times10^{-17}$',r'$z>1.0, [OII]>8\times10^{-17}$'],
    					     ['magenta','tan','powderblue','blue'],
    					     ['^','s','o','o']):
    	    cut= keep[cut_name]
    	    ax.scatter(tab.rz[cut],tab.gr[cut], 
    		       marker=marker, color=color, label=lab)
    	    ax.set_xlim(xrange)
    	    ax.set_ylim(yrange)
    	    xlab= ax.set_xlabel('r-z')
    	    ylab= ax.set_ylabel('g-r')
    	    leg=ax.legend(loc=(0,1.05), ncol=2,prop={'size': 14}, labelspacing=0.2,
    			  markerscale=1.5)
    	    name='dr%d_FDR_ELG.png' % self.DR
    	    kwargs= dict(bbox_extra_artists=[leg,xlab,ylab], bbox_inches='tight',dpi=150)
    	if self.savefig:
    	    plt.savefig(name, **kwargs)
    	    plt.close()
    	    print('Wrote {}'.format(name))

    def plot_FDR_multi(self):
    	tab= self.get_dr3_deep2()
    	# Cuts 'lowz','medz_lowO2' ...
    	keep= self.get_FDR_cuts(tab)
    	# Plot
    	fig,ax = plt.subplots(1,4,sharex=True,sharey=True,figsize=(18,4))
    	plt.subplots_adjust(wspace=0.1,hspace=0)
    	ts= TSBox(src='ELG')
    	xrange,yrange= xyrange['x_elg'],xyrange['y_elg']
    	for cnt,cut_name,lab,color,marker in zip(range(4),
    						 ['lowz','medz_lowO2','medz_hiO2','hiz_hiO2'],
    						 [r'$z<0.6$',r'$z>0.6, [OII]<8\times10^{-17}$',
    						  r'$z>0.6, [OII]>8\times10^{-17}$',r'$z>1.0, [OII]>8\times10^{-17}$'],
    						 ['magenta','tan','powderblue','blue'],
    						 ['^','s','o','o']):
                # Add box
    	    ts.add_ts_box(ax[cnt], xlim=xrange,ylim=yrange)
    	    # Add points
    	    cut= keep[cut_name]
    	    ax[cnt].scatter(tab.rz[cut],tab.gr[cut], marker=marker,color=color)
    	    ti_loc=ax[cnt].set_title(lab)
    	    ax[cnt].set_xlim(xrange)
    	    ax[cnt].set_ylim(yrange)
    	    xlab= ax[cnt].set_xlabel('r-z')
    	    ylab= ax[cnt].set_ylabel('g-r')
    	name='dr%d_ELG_FDR_multi.png' % self.DR
    	kwargs= dict(bbox_extra_artists=[ti_loc,xlab,ylab], bbox_inches='tight',dpi=150)
    	if self.savefig:
    	    plt.savefig(name, **kwargs)
    	    plt.close()
    	    print('Wrote {}'.format(name))
        
    def plot_obiwan_multi(self):
    	tab= self.get_dr3_deep2()
    	# Cuts 'lowz','medz_lowO2' ...
    	keep= self.get_obiwan_cuts(tab)
    	# Plot
    	fig,ax = plt.subplots(1,2,sharex=True,sharey=True,figsize=(9,4))
    	plt.subplots_adjust(wspace=0.1,hspace=0)
    	ts= TSBox(src='ELG')
    	xrange,yrange= xyrange['x_elg'],xyrange['y_elg']
    	for cnt,thecut,lab,color in zip(range(2),
    					[keep, keep == False],
    					[r'$0.8<z<1.4, [OII] > 0$','Everything else'],
    					['b','g']):
    	    # Add box
    	    ts.add_ts_box(ax[cnt], xlim=xrange,ylim=yrange)
    	    # Add points
    	    ax[cnt].scatter(tab.rz[thecut],tab.gr[thecut], marker='o',color=color)
    	    ti_loc=ax[cnt].set_title(lab)
    	    ax[cnt].set_xlim(xrange)
    	    ax[cnt].set_ylim(yrange)
    	    xlab= ax[cnt].set_xlabel('r-z')
    	    ylab= ax[cnt].set_ylabel('g-r')
    	    name='dr%d_ELG_obiwan_multi.png' % self.DR
    	    kwargs= dict(bbox_extra_artists=[ti_loc,xlab,ylab], bbox_inches='tight',dpi=150)
    	    if self.savefig:
    		plt.savefig(name, **kwargs)
    		plt.close()
    		print('Wrote {}'.format(name))
    
    def plot_LRG_FDR_wELG_data(self):
    	dic= self.get_elgs_FDR_cuts()
    	ts= TSBox(src='LRG')
    	xrange,yrange= xyrange['x_lrg'],xyrange['y_lrg']
    	# Plot
    	fig,ax = plt.subplots(1,4,sharex=True,sharey=True,figsize=(18,4))
    	plt.subplots_adjust(wspace=0.1,hspace=0)
    	for cnt,key,col,marker,ti in zip(range(4),\
    					 ['loz','oiifaint','oiibright_loz','oiibright_hiz'],\
    					 ['magenta','tan','powderblue','blue'],\
    					 ['^','s','o','o'],\
    					 [r'$z<0.6$',r'$z>0.6, [OII]<8\times10^{-17}$',r'$z>0.6, [OII]>8\times10^{-17}$',r'$z>1.0, [OII]>8\times10^{-17}$']):
    	    # Add box
    	    ts.add_ts_box(ax[cnt], xlim=xrange,ylim=yrange)
    	    # Add points
    	    b= key
    	    ax[cnt].scatter(dic['rz'][b],dic['rw1'][b], marker=marker,color=col)
    	    ti_loc=ax[cnt].set_title(ti)
    	    ax[cnt].set_xlim(xrange)
    	    ax[cnt].set_ylim(yrange)
    	    xlab= ax[cnt].set_xlabel('r-z')
    	    ylab= ax[cnt].set_ylabel('r-W1')
    	name='dr%d_LRG_FDR_wELG_data.png' % self.DR
    	kwargs= dict(bbox_extra_artists=[ti_loc,xlab,ylab], bbox_inches='tight',dpi=150)
    	if self.savefig:
    	    plt.savefig(name, **kwargs)
    	    plt.close()
    	    print('Wrote {}'.format(name))
            
	def plot_FDR_mag_dist(self):
		tab= self.get_dr3_deep2()
		keep= self.get_FDR_cuts(tab)
		# Plot
		fig,ax = plt.subplots(1,4,sharex=True,sharey=True,figsize=(18,4))
		plt.subplots_adjust(wspace=0.1,hspace=0)
		for cnt,cut_name,lab in zip(range(4),
						['lowz','medz_lowO2','medz_hiO2','hiz_hiO2'],
						[r'$z<0.6$',r'$z>0.6, [OII]<8\times10^{-17}$',
							r'$z>0.6, [OII]>8\times10^{-17}$',r'$z>1.0, [OII]>8\times10^{-17}$']):
			# Mag data
			cut= keep[cut_name]
			mags=dict(r= tab.r_wdust[cut],
					  g= tab.gr[cut] + tab.r_wdust[cut],
					  z= tab.r_wdust[cut] - tab.rz[cut]) 
			for band,color in zip(['g','r','z'],['g','r','m']):
				# nans present
				mags[band]= mags[band][ np.isfinite(mags[band]) ]
				# histograms
				h,edges= np.histogram(mags[band],bins=40,normed=True)
				binc= (edges[1:]+edges[:-1])/2.
				ax[cnt].step(binc,h,where='mid',lw=1,c=color,label='%s' % band)
			ti_loc=ax[cnt].set_title(lab)
			ax[cnt].set_xlim([20,26])
			ax[cnt].set_ylim([0,0.9])
			xlab= ax[cnt].set_xlabel('AB mag')
			ylab= ax[cnt].set_ylabel('PDF')
			name='dr%d_ELG_mag_dist_FDR.png' % self.DR
		kwargs= dict(bbox_extra_artists=[ti_loc,xlab,ylab], bbox_inches='tight',dpi=150)
		if self.savefig:
			plt.savefig(name, **kwargs)
			plt.close()
			print('Wrote {}'.format(name))

	def plot_obiwan_mag_dist(self):
		tab= self.get_dr3_deep2()
		keep= self.get_obiwan_cuts(tab)
		# Plot
		fig,ax = plt.subplots(1,2,sharex=True,sharey=True,figsize=(9,4))
		plt.subplots_adjust(wspace=0.1,hspace=0)
		for cnt,thecut,lab in zip(range(4),
						[keep, keep == False],
						[r'$0.8<z<1.4, [OII] > 0$','Everything else']):
			# Mag data
			mags=dict(r= tab.r_wdust[thecut],
					  g= tab.gr[thecut] + tab.r_wdust[thecut],
					  z= tab.r_wdust[thecut] - tab.rz[thecut]) 
			for band,color in zip(['g','r','z'],['g','r','m']):
				# nans present
				mags[band]= mags[band][ np.isfinite(mags[band]) ]
				# histograms
				h,edges= np.histogram(mags[band],bins=40,normed=True)
				binc= (edges[1:]+edges[:-1])/2.
				ax[cnt].step(binc,h,where='mid',lw=1,c=color,label='%s' % band)
			ti_loc=ax[cnt].set_title(lab)
			ax[cnt].set_xlim([20,26])
			ax[cnt].set_ylim([0,0.9])
			xlab= ax[cnt].set_xlabel('AB mag')
			ylab= ax[cnt].set_ylabel('PDF')
			name='dr%d_ELG_mag_dist_obiwan.png' % self.DR
		kwargs= dict(bbox_extra_artists=[ti_loc,xlab,ylab], bbox_inches='tight',dpi=150)
		if self.savefig:
			plt.savefig(name, **kwargs)
			plt.close()
			print('Wrote {}'.format(name))


	def plot(self):
		self.plot_FDR()
		self.plot_FDR_multipanel()
		#plot_FDR(self.Xall,self.cuts,src='ELG')
		#b= self.cuts['any_elg']
		#color_color_plot(self.Xall[b,:],src='ELG',append='_FDR') #,extra=True)
		#Xall,cuts, morph= elg_data()
		#color_color_plot(Xall, src='ELG',append='_synth') #,extra=True)
		#b= cuts['has_morph']
		#color_color_plot(Xall[b,:],src='ELG',append='_synth+morph') #,extra=True)

	def plot_kde(self):
		rz,gr,r_nodust,r_wdust,redshift= self.get_elgs_FDR_cuts()
		x= r_wdust['med2hiz_oiibright']
		y= rz['med2hiz_oiibright']
		z= gr['med2hiz_oiibright']
		d4= redshift['med2hiz_oiibright']
		cut= (np.isfinite(x))* (np.isfinite(y))* (np.isfinite(z))
		x,y,z,d4= x[cut],y[cut],z[cut],d4[cut]
		labels=['r wdust','r-z','g-r','redshift']
		kde_obj= KernelOfTruth([x,y,z,d4],labels,\
							   [(20.5,25.),(0,2),(-0.5,1.5),(0.6,1.6)],\
							   bandwidth=0.05,
							   kdefn=self.kdefn,loadkde=self.loadkde)
		xylims=dict(x1=(20.5,25.5),y1=(0,0.8),\
					x2=xyrange['x_elg'],y2=xyrange['y_elg'],\
					x3=(0.6,1.6),y3=(0.,1.0))
		#kde_obj.plot_1band_and_color(ndraws=1000,xylims=xylims,prefix='elg_')
		kde_obj.plot_1band_color_and_redshift(ndraws=1000,xylims=xylims,prefix='elg_')
		if self.savekde:
			if os.path.exists(self.kdefn):
				os.remove(self.kdefn)
			kde_obj.save(name=self.kdefn)

	def plot_kde_shapes(self):
		re,n,ba,pa= self.get_acs_matched_deep2()
		pa+= 90. # 0-180 deg
		#cut= (np.isfinite(x))* (np.isfinite(y))* (np.isfinite(z))
		#x,y,z,d4= x[cut],y[cut],z[cut],d4[cut]
		# ba > 0
		labels=['re','n','ba','pa']
		kde_obj= KernelOfTruth([re,n,ba,pa],labels,\
							   [(0.,100.),(0.,10.),(0.2,0.9),(0.,180.)],\
							   bandwidth=0.05,kernel='tophat',\
							   kdefn=self.kde_shapes_fn,loadkde=self.loadkde)
		xylims=dict(x1=(0,100),\
					x2=(0,10),\
					x3=(0,1),\
					x4=(0,180))
		#kde_obj.plot_1band_and_color(ndraws=1000,xylims=xylims,prefix='elg_')
		kde_obj.plot_galaxy_shapes(ndraws=1000,xylims=xylims,name='elg_shapes_kde.png')
		if self.savekde:
			if os.path.exists(self.kde_shapes_fn):
				os.remove(self.kde_shapes_fn)
			kde_obj.save(name=self.kde_shapes_fn)


class FDR_lrg(object):
	def get_FDR_cuts(self,tab):
		keep={}  
		keep['star']= tab.type == 1
		keep['blue_galaxy']= (tab.type == 0) *\
							 (tab.mod_gal > 8)
		keep['red_galaxy_lowz']= (tab.type == 0) * \
								 (tab.mod_gal <= 8) *\
								 (tab.zp_gal <= 0.6)
		keep['red_galaxy_hiz']= (tab.type == 0) * \
								 (tab.mod_gal <= 8) *\
								 (tab.zp_gal > 0.6)
		return keep

                        
	def plot_FDR(self):
		tab= self.get_dr3_cosmos()
		keep= self.get_FDR_cuts(tab)
		# Plot
		fig,ax = plt.subplots()
		rgb_cols=get_rgb_cols()
		# Add box
		ts= TSBox(src='LRG')
		xrange,yrange= xyrange['x_lrg'],xyrange['y_lrg']
		ts.add_ts_box(ax, xlim=xrange,ylim=yrange)
		# Data
		for cnt,cut_name,rgb in zip(range(4),
						['star','red_galaxy_lowz','red_galaxy_hiz','blue_galaxy'],
									rgb_cols):
			rgb= (rgb[0]/255.,rgb[1]/255.,rgb[2]/255.)
			cut= keep[cut_name]
			ax.scatter(tab.rz[cut],tab.rW1[cut],c=[rgb],
					   edgecolors='none',marker='o',s=10.,rasterized=True,label=cut_name)
		ax.set_xlim(xrange)
		ax.set_ylim(yrange)
		xlab=ax.set_xlabel('r-z')
		ylab=ax.set_ylabel('r-W1')
		leg=ax.legend(loc=(0,1.05), ncol=2,prop={'size': 14}, labelspacing=0.2,\
					  markerscale=2,scatterpoints=1)
		#handles,labels = ax.get_legend_handles_labels()
		#index=[0,1,2,3]
		#handles,labels= np.array(handles)[index],np.array(labels)[index]
		#leg=ax.legend(handles,labels,loc=(0,1.05),ncol=2,scatterpoints=1,markerscale=2)
		name='dr%d_FDR_LRG.png' % self.DR
		if self.savefig:
			plt.savefig(name,\
						bbox_extra_artists=[leg,xlab,ylab], bbox_inches='tight',dpi=150)
			plt.close()
			print('Wrote {}'.format(name))

	def plot_FDR_multi(self):
		tab= self.get_dr3_cosmos()
		keep= self.get_FDR_cuts(tab)
		# Plot
		fig,ax = plt.subplots(1,4,sharex=True,sharey=True,figsize=(16,4))
		plt.subplots_adjust(wspace=0.1,hspace=0)
		rgb_cols=get_rgb_cols()
		for cnt,cut_name,rgb in zip(range(4),   
					   ['star','red_galaxy_lowz','red_galaxy_hiz','blue_galaxy'],
									rgb_cols):
			rgb= (rgb[0]/255.,rgb[1]/255.,rgb[2]/255.)
			cut= keep[cut_name]
			ax[cnt].scatter(tab.rz[cut],tab.rW1[cut],c=[rgb],
							edgecolors='none',marker='o',s=10.,rasterized=True)#,label=key)
			ti=ax[cnt].set_title(cut_name)
			ax[cnt].set_xlim([0,2.5])
			ax[cnt].set_ylim([-2,6])
			xlab=ax[cnt].set_xlabel('r-z')
			# Add box
			ts= TSBox(src='LRG')
			xrange,yrange= xyrange['x_lrg'],xyrange['y_lrg']
			ts.add_ts_box(ax[cnt], xlim=xrange,ylim=yrange)
			ylab=ax[cnt].set_ylabel('r-W1')
		#handles,labels = ax.get_legend_handles_labels()
		#index=[0,1,2,3]
		#handles,labels= np.array(handles)[index],np.array(labels)[index]
		#leg=ax.legend(handles,labels,loc=(0,1.05),ncol=2,scatterpoints=1,markerscale=2)
		name='dr%d_LRG_FDR_multi.png' % self.DR
		if self.savefig:
			plt.savefig(name,\
						bbox_extra_artists=[ti,xlab,ylab], bbox_inches='tight',dpi=150)
			plt.close()
			print('Wrote {}'.format(name))

	def plot_obiwan_multi(self):
		tab= self.get_dr3_cosmos()
		keep= self.get_obiwan_cuts(tab)
		# Plot
		fig,ax = plt.subplots(1,2,sharex=True,sharey=True,figsize=(8,4))
		plt.subplots_adjust(wspace=0.1,hspace=0)
		rgb_cols=get_rgb_cols()
		for cnt,thecut,cut_name,rgb in zip(range(2),   
								  [keep, keep == False],
								  ['red galaxy','everything else'],
								  rgb_cols):
			rgb= (rgb[0]/255.,rgb[1]/255.,rgb[2]/255.)
			ax[cnt].scatter(tab.rz[thecut],tab.rW1[thecut],c=[rgb],
							edgecolors='none',marker='o',s=10.,rasterized=True)#,label=key)
			ti=ax[cnt].set_title(cut_name)
			ax[cnt].set_xlim([0,2.5])
			ax[cnt].set_ylim([-2,6])
			xlab=ax[cnt].set_xlabel('r-z')
			# Add box
			ts= TSBox(src='LRG')
			xrange,yrange= xyrange['x_lrg'],xyrange['y_lrg']
			ts.add_ts_box(ax[cnt], xlim=xrange,ylim=yrange)
			ylab=ax[cnt].set_ylabel('r-W1')
		#handles,labels = ax.get_legend_handles_labels()
		#index=[0,1,2,3]
		#handles,labels= np.array(handles)[index],np.array(labels)[index]
		#leg=ax.legend(handles,labels,loc=(0,1.05),ncol=2,scatterpoints=1,markerscale=2)
		name='dr%d_LRG_obiwan_multi.png' % self.DR
		if self.savefig:
			plt.savefig(name,\
						bbox_extra_artists=[ti,xlab,ylab], bbox_inches='tight',dpi=150)
			plt.close()
			print('Wrote {}'.format(name))


	def plot_LRGs_in_ELG_FDR(self):
		data= self.get_lrgs_FDR_cuts()
		fig,ax = plt.subplots(1,4,sharex=True,sharey=True,figsize=(16,4))
		plt.subplots_adjust(wspace=0.1,hspace=0)
		rgb_cols=get_rgb_cols()
		keys= ['star','red_galaxy_lowz','red_galaxy_hiz','blue_galaxy']
		for cnt,key,rgb in zip(range(4),keys,rgb_cols):
			rgb= (rgb[0]/255.,rgb[1]/255.,rgb[2]/255.)
			ax[cnt].scatter(data['rz'][key],data['g_wdust'][key]-data['r_wdust'][key],c=[rgb],
							edgecolors='none',marker='o',s=10.,rasterized=True)#,label=key)
			ti=ax[cnt].set_title(key)
			xrange,yrange= xyrange['x_elg'],xyrange['y_elg']
			ax[cnt].set_xlim(xrange)
			ax[cnt].set_ylim(yrange)
			xlab=ax[cnt].set_xlabel('r-z')
			ylab=ax[cnt].set_ylabel('g-r')
			# Add box
			ts= TSBox(src='ELG')
			ts.add_ts_box(ax[cnt], xlim=xrange,ylim=yrange)
		#handles,labels = ax.get_legend_handles_labels()
		#index=[0,1,2,3]
		#handles,labels= np.array(handles)[index],np.array(labels)[index]
		#leg=ax.legend(handles,labels,loc=(0,1.05),ncol=2,scatterpoints=1,markerscale=2)
		name='dr%d_LRGs_in_ELG_FDR.png' % self.DR
		if self.savefig:
			plt.savefig(name,\
						bbox_extra_artists=[ti,xlab,ylab], bbox_inches='tight',dpi=150)
			plt.close()
			print('Wrote {}'.format(name))

	def plot_FDR_mag_dist(self):
		tab= self.get_dr3_cosmos()
		keep= self.get_FDR_cuts(tab)
		fig,ax = plt.subplots(1,4,sharex=True,sharey=True,figsize=(16,4))
		plt.subplots_adjust(wspace=0.1,hspace=0)
		rgb_cols=get_rgb_cols()
		for cnt,cut_name,rgb in zip(range(4),
					['star','red_galaxy_lowz','red_galaxy_hiz','blue_galaxy'],
									rgb_cols):
			rgb= (rgb[0]/255.,rgb[1]/255.,rgb[2]/255.)
			# Mag data
			cut= keep[cut_name]
			mags=dict(r= tab.r_wdust[cut],
					  g= tab.g_wdust[cut],
					  z= tab.z_wdust[cut]) 
			for band,color in zip(['g','r','z'],['g','r','m']):
				# nans present
				mags[band]= mags[band][ np.isfinite(mags[band]) ]
				# histograms
				h,edges= np.histogram(mags[band],bins=40,normed=True)
				binc= (edges[1:]+edges[:-1])/2.
				ax[cnt].step(binc,h,where='mid',lw=1,c=color,label='%s' % band)
			ti=ax[cnt].set_title(cut_name)
			ax[cnt].set_xlim([16,26])
			ax[cnt].set_ylim([0,0.9])
			xlab=ax[cnt].set_xlabel('AB mags')
			ylab=ax[cnt].set_ylabel('PDF')
		#handles,labels = ax.get_legend_handles_labels()
		#index=[0,1,2,3]
		#handles,labels= np.array(handles)[index],np.array(labels)[index]
		#leg=ax.legend(handles,labels,loc=(0,1.05),ncol=2,scatterpoints=1,markerscale=2)
		name='dr%d_LRG_mag_dist_FDR.png' % self.DR
		if self.savefig:
			plt.savefig(name,\
						bbox_extra_artists=[ti,xlab,ylab], bbox_inches='tight',dpi=150)
			plt.close()
			print('Wrote {}'.format(name))

	def plot_obiwan_mag_dist(self):
		tab= self.get_dr3_cosmos()
		keep= self.get_obiwan_cuts(tab)
		fig,ax = plt.subplots(1,2,sharex=True,sharey=True,figsize=(8,4))
		plt.subplots_adjust(wspace=0.1,hspace=0)
		rgb_cols=get_rgb_cols()
		for cnt,thecut,lab,rgb in zip(range(2),
					[keep, keep == False],
					['red galaxy','everything else'],
									rgb_cols):
			rgb= (rgb[0]/255.,rgb[1]/255.,rgb[2]/255.)
			# Mag data
			mags=dict(r= tab.r_wdust[thecut],
					  g= tab.g_wdust[thecut],
					  z= tab.z_wdust[thecut]) 
			for band,color in zip(['g','r','z'],['g','r','m']):
				# nans present
				mags[band]= mags[band][ np.isfinite(mags[band]) ]
				# histograms
				h,edges= np.histogram(mags[band],bins=40,normed=True)
				binc= (edges[1:]+edges[:-1])/2.
				ax[cnt].step(binc,h,where='mid',lw=1,c=color,label='%s' % band)
			ti=ax[cnt].set_title(lab)
			ax[cnt].set_xlim([16,26])
			ax[cnt].set_ylim([0,0.9])
			xlab=ax[cnt].set_xlabel('AB mags')
			ylab=ax[cnt].set_ylabel('PDF')
		#handles,labels = ax.get_legend_handles_labels()
		#index=[0,1,2,3]
		#handles,labels= np.array(handles)[index],np.array(labels)[index]
		#leg=ax.legend(handles,labels,loc=(0,1.05),ncol=2,scatterpoints=1,markerscale=2)
		name='dr%d_LRG_mag_dist_obiwan.png' % self.DR
		if self.savefig:
			plt.savefig(name,\
						bbox_extra_artists=[ti,xlab,ylab], bbox_inches='tight',dpi=150)
			plt.close()
			print('Wrote {}'.format(name))

	def plot_vipers(self):
		rz,rW1= self.get_vipers()
		fig,ax = plt.subplots(figsize=(5,4))
		plt.subplots_adjust(wspace=0,hspace=0)
		rgb=get_rgb_cols()[0]
		rgb= (rgb[0]/255.,rgb[1]/255.,rgb[2]/255.)
		ax.scatter(rz,rW1,c=[rgb],edgecolors='none',marker='o',s=10.,rasterized=True)#,label=key)
		ax.set_xlim([0,2.5])
		ax.set_ylim([-2,6])
		xlab=ax.set_xlabel('r-z')
		ylab=ax.set_ylabel('r-W1')
		# Add box
		ts= TSBox(src='LRG')
		xrange,yrange= xyrange['x_lrg'],xyrange['y_lrg']
		ts.add_ts_box(ax, xlim=xrange,ylim=yrange)
		#handles,labels = ax.get_legend_handles_labels()
		#index=[0,1,2,3]
		#handles,labels= np.array(handles)[index],np.array(labels)[index]
		#leg=ax.legend(handles,labels,loc=(0,1.05),ncol=2,scatterpoints=1,markerscale=2)
		name='dr%d_LRG_vipers.png' % self.DR
		if self.savefig:
			plt.savefig(name,\
						bbox_extra_artists=[xlab,ylab], bbox_inches='tight',dpi=150)
			plt.close()
			print('Wrote {}'.format(name))


	def plot(self):
		self.plot_FDR()
		self.plot_FDR_multipanel()
		self.plot_vipers()
		#plot_FDR(self.Xall,self.cuts,src='LRG')
		#color_color_plot(self.Xall,src='LRG',append='cc') #,extra=True)
		#b= self.cuts['lrg')
		#color_color_plot(self.Xall[b,:],src='LRG') #,extra=True)

	def plot_kde(self,loadkde=False,savekde=False):
		"""No Targeting cuts on g band, but need to fit it so can insert in grz image"""
		rz,rW1,r_nodust,r_wdust,z_nodust,z_wdust,g_wdust,redshift= self.get_lrgs_FDR_cuts()
		x= z_wdust['red_galaxy']
		y= rz['red_galaxy']
		z= rW1['red_galaxy']
		d4= redshift['red_galaxy']
		M= g_wdust['red_galaxy']
		cut= (np.isfinite(x))* (np.isfinite(y))* (np.isfinite(z))* (np.isfinite(M))
		# Redshift > 0 given bandwidth
		bandwidth=0.05
		cut*= (d4 - bandwidth >= 0.)
		x,y,z,d4,M= x[cut],y[cut],z[cut],d4[cut],M[cut]
		labels=['z wdust','r-z','r-W1','redshift','g wdust']
		kde_obj= KernelOfTruth([x,y,z,d4,M],labels,\
						   [(17.,22.),(0,2.5),(-2,5.),(0.,1.6),(17.,29)],\
						   bandwidth=bandwidth,kernel='tophat',\
						   kdefn=self.kdefn,loadkde=self.loadkde)
		xylims=dict(x1=(17.,22.),y1=(0,0.7),\
					x2=xyrange['x_lrg'],y2=xyrange['y_lrg'],\
					x3=(0.,1.6),y3=(0,1.),\
					x4=(17.,29),y4=(0,0.7))
		#kde_obj.plot_1band_and_color(ndraws=1000,xylims=xylims,prefix='lrg_')
		kde_obj.plot_1band_color_and_redshift(ndraws=1000,xylims=xylims,prefix='lrg_')
		if self.savekde:
			if os.path.exists(self.kdefn):
				os.remove(self.kdefn)
			kde_obj.save(name=self.kdefn)


	def plot_kde_shapes(self):
		re,n,ba,pa= self.get_acs_matched_cosmoszphot()
		pa+= 90. # 0-180 deg
		#cut= (np.isfinite(x))* (np.isfinite(y))* (np.isfinite(z))
		#x,y,z,d4= x[cut],y[cut],z[cut],d4[cut]
		# ba > 0
		labels=['re','n','ba','pa']
		kde_obj= KernelOfTruth([re,n,ba,pa],labels,\
							   [(0.,100.),(0.,10.),(0.2,0.9),(0.,180.)],\
							   bandwidth=0.05,kernel='tophat',\
							   kdefn=self.kde_shapes_fn,loadkde=self.loadkde)
		xylims=dict(x1=(-10,100),\
					x2=(-2,10),\
					x3=(-0.2,1.2),\
					x4=(-20,200))
		#kde_obj.plot_1band_and_color(ndraws=1000,xylims=xylims,prefix='elg_')
		kde_obj.plot_galaxy_shapes(ndraws=10000,xylims=xylims,name='lrg_shapes_kde.png')
		if self.savekde:
			if os.path.exists(self.kde_shapes_fn):
				os.remove(self.kde_shapes_fn)
			kde_obj.save(name=self.kde_shapes_fn)

class FDR_qso(object):
	def plot_FDR(self):
		# Data
		qsos= self.get_qsos()
		star_obj= STAR(DR=self.DR,savefig=False)
		stars= star_obj.get_purestars()
		hiz=2.1
		index={}
		index['hiz']= qsos.get('z') > hiz
		index['loz']= qsos.get('z') <= hiz
		# Plot
		fig,ax = plt.subplots(1,2,figsize=(10,4))
		plt.subplots_adjust(wspace=0.1,hspace=0)
		# Stars
		ax[0].scatter(stars.get('decam_mag_nodust')[:,2]-stars.get('decam_mag_nodust')[:,4],\
					  stars.get('decam_mag_nodust')[:,1]-stars.get('decam_mag_nodust')[:,2],\
					  c='b',edgecolors='none',marker='o',s=10.,rasterized=True, label='stars',alpha=self.alpha)
		W= 0.75*stars.get('wise_mag_nodust')[:,0]+ 0.25*stars.get('wise_mag_nodust')[:,1]
		ax[1].scatter(stars.get('decam_mag_nodust')[:,1]-stars.get('decam_mag_nodust')[:,4],\
					  stars.get('decam_mag_nodust')[:,2]-W,\
					  c='b',edgecolors='none',marker='o',s=10.,rasterized=True, label='stars',alpha=self.alpha)
		# QSOs
		for key,lab,col in zip(['loz','hiz'],['(z < 2.1)','(z > 2.1)'],['magenta','red']):
			i= index[key]
			ax[0].scatter(qsos.get('decam_mag_nodust')[:,2][i]-qsos.get('decam_mag_nodust')[:,4][i],\
						  qsos.get('decam_mag_nodust')[:,1][i]-qsos.get('decam_mag_nodust')[:,2][i],\
						  c=col,edgecolors='none',marker='o',s=10.,rasterized=True, label='qso '+lab,alpha=self.alpha)
			W= 0.75*qsos.get('wise_mag_nodust')[:,0]+ 0.25*qsos.get('wise_mag_nodust')[:,1]
			ax[1].scatter(qsos.get('decam_mag_nodust')[:,1][i]-qsos.get('decam_mag_nodust')[:,4][i],\
						  qsos.get('decam_mag_nodust')[:,2][i]-W[i],\
						  c=col,edgecolors='none',marker='o',s=10.,rasterized=True, label='qso '+lab,alpha=self.alpha)
		
		#for xlim,ylim,x_lab,y_lab in ax[0].set_xlim([-0.5,3.])
		ax[0].set_xlim(xyrange['x1_qso'])
		ax[1].set_xlim(xyrange['x2_qso'])
		ax[0].set_ylim(xyrange['y1_qso'])
		ax[1].set_ylim(xyrange['y2_qso'])
		xlab=ax[0].set_xlabel('r-z')
		xlab=ax[1].set_xlabel('g-z')
		ylab=ax[0].set_ylabel('g-r')
		ylab=ax[1].set_ylabel('r-W')
		leg=ax[0].legend(loc=(0,1.02),scatterpoints=1,ncol=3,markerscale=2)
		## Add box
		#ts= TSBox(src='LRG')
		#xrange,yrange= xyrange['x_lrg'],xyrange['y_lrg']
		#ts.add_ts_box(ax[cnt], xlim=xrange,ylim=yrange)
		name='dr%d_FDR_QSO.png' % self.DR
		if self.savefig:
			plt.savefig(name,\
						bbox_extra_artists=[leg,xlab,ylab], bbox_inches='tight',dpi=150)
			plt.close()
			print('Wrote {}'.format(name))

	def plot_FDR_multipanel(self):
		# Data
		qsos= self.get_qsos()
		star_obj= STAR(DR=self.DR,savefig=False)
		stars= star_obj.get_purestars()
		hiz=2.1
		index={}
		index['hiz']= qsos.get('z') > hiz
		index['loz']= qsos.get('z') <= hiz
		# Plot
		fig,ax = plt.subplots(3,2,figsize=(10,12))
		plt.subplots_adjust(wspace=0.2,hspace=0.1)
		# Stars top panel
		ax[0,0].scatter(stars.get('decam_mag_nodust')[:,2]-stars.get('decam_mag_nodust')[:,4],\
					  stars.get('decam_mag_nodust')[:,1]-stars.get('decam_mag_nodust')[:,2],\
					  c='b',edgecolors='none',marker='o',s=10.,rasterized=True, label='stars',alpha=self.alpha)
		W= 0.75*stars.get('wise_mag_nodust')[:,0]+ 0.25*stars.get('wise_mag_nodust')[:,1]
		ax[0,1].scatter(stars.get('decam_mag_nodust')[:,1]-stars.get('decam_mag_nodust')[:,4],\
					  stars.get('decam_mag_nodust')[:,2]-W,\
					  c='b',edgecolors='none',marker='o',s=10.,rasterized=True, label='stars',alpha=self.alpha)
		# QSOs loz middle, hiz bottom
		for cnt,key,lab,col in zip([1,2],['loz','hiz'],['(z < 2.1)','(z > 2.1)'],['magenta','red']):
			i= index[key]
			ax[cnt,0].scatter(qsos.get('decam_mag_nodust')[:,2][i]-qsos.get('decam_mag_nodust')[:,4][i],\
						  qsos.get('decam_mag_nodust')[:,1][i]-qsos.get('decam_mag_nodust')[:,2][i],\
						  c=col,edgecolors='none',marker='o',s=10.,rasterized=True, label='qso '+lab,alpha=self.alpha)
			W= 0.75*qsos.get('wise_mag_nodust')[:,0]+ 0.25*qsos.get('wise_mag_nodust')[:,1]
			ax[cnt,1].scatter(qsos.get('decam_mag_nodust')[:,1][i]-qsos.get('decam_mag_nodust')[:,4][i],\
						  qsos.get('decam_mag_nodust')[:,2][i]-W[i],\
						  c=col,edgecolors='none',marker='o',s=10.,rasterized=True, label='qso '+lab,alpha=self.alpha)
		
		for cnt in range(3):
			#for xlim,ylim,x_lab,y_lab in ax[0].set_xlim([-0.5,3.])
			ax[cnt,0].set_xlim(xyrange['x1_qso'])
			ax[cnt,1].set_xlim(xyrange['x2_qso'])
			ax[cnt,0].set_ylim(xyrange['y1_qso'])
			ax[cnt,1].set_ylim(xyrange['y2_qso'])
			xlab=ax[cnt,0].set_xlabel('r-z')
			xlab=ax[cnt,1].set_xlabel('g-z')
			ylab=ax[cnt,0].set_ylabel('g-r')
			ylab=ax[cnt,1].set_ylabel('r-W')
		#leg=ax[0,0].legend(loc=(0,1.02),scatterpoints=1,ncol=3,markerscale=2)
		## Add box
		#ts= TSBox(src='LRG')
		#xrange,yrange= xyrange['x_lrg'],xyrange['y_lrg']
		#ts.add_ts_box(ax[cnt], xlim=xrange,ylim=yrange)
		name='dr%d_FDR_QSO_multi.png' % self.DR
		if self.savefig:
			plt.savefig(name,\
						bbox_extra_artists=[leg,xlab,ylab], bbox_inches='tight',dpi=150)
			plt.close()
			print('Wrote {}'.format(name))

	def plot(self):
		self.plot_FDR()
		self.plot_FDR_multipanel()

                        
                        


def plot_tractor_galfit_shapes(cat,prefix=''):
    for name in ['re','n']:
        fig,ax= plt.subplots()
        keep= cat.get('tractor_'+name) > 0.
        ax.scatter(cat.get(name)[keep], cat.get('tractor_'+name)[keep],
                   c='b',marker='o',s=10.,rasterized=True)
        xlab= ax.set_xlabel('galfit_hi_%s' % name)
        ylab= ax.set_ylabel('tractor_%s' % name)
        savenm= 'tractor_galfit_%s_%s.png' % (prefix,name)
        plt.savefig(savenm,bbox_extra_artists=[xlab,ylab], bbox_inches='tight',dpi=150)
        plt.close()
        print('Wrote %s' % savenm)
 
# 2D plots




class GalaxyPrior(object):
	"""Stores all 4 galaxy type classes"""

	def __init__(self):
		#self.qso= QSO()
		self.star= STAR()
		self.lrg= LRG()
		self.elg= ELG()
	def plot_all(self):
		#self.qso.plot()
		self.star.plot()
		self.lrg.plot()
		self.elg.plot()

def get_parser():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,\
                                     description='Generate a legacypipe-compatible CCDs file \
                                                  from a set of reduced imaging.')
    parser.add_argument('--outdir',action='store',default='./',help='where to write KDE fits and QA plots',required=False)
    return parser
 

if __name__ == '__main__':
    import argparse
    parser= get_parser()
    args = parser.parse_args()

    elg= EmptyClass()
    
    d= Data()
    outdir='/home/kaylan/mydata/priors_data'
    d.fetch(outdir)
    elg.data= d.load_elg(DR=3)

    elg.model= KDE_Model('elg',elg.data,outdir)
    elg.model.kde= elg.model.get_kde()

    p= Plot(outdir)
    p.elg(elg.data, elg.model.df_wcut,elg.model.df_for_kde,
          elg.model.kde)
    raise ValueError('DONE')

    #elg.plot_FDR()
    #elg.plot_FDR_multi()
    #elg.plot_obiwan_multi()
    #elg.plot_FDR_mag_dist()
    #elg.plot_obiwan_mag_dist()
    #elg.plot_LRG_FDR_wELG_data()
    #elg.plot_FDR_multipanel()
    #elg.get_acs_matched_deep2()
    #elg.plot_dr3_acs_deep2()
    #elg.plot_kde_shapes()
    #elg.plot_kde()
    #elg.plot_redshift()
    #elg.cross_validate_redshift()

    #elg.plot_kde()
    #elg.plot()
    #lrg.plot_dr3_cosmos_acs()
    #lrg.plot_FDR()
    #lrg.plot_FDR_multi()
    #lrg.plot_obiwan_multi()
    #lrg.plot_FDR_mag_dist()
    #lrg.plot_obiwan_mag_dist()
    #lrg.plot_LRG_FDR_mag_dist()
    #lrg.plot_LRGs_in_ELG_FDR()
    #lrg.plot_FDR()
    #lrg.plot_FDR_multipanel()
    #lrg.plot_kde_shapes()
    #lrg.plot_kde()
    #lrg.plot()

    

    

    
