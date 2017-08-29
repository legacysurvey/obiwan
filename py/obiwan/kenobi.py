# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""This is obiwan.kenobi
...it runs the legacypipe/Tractor pipeline on images containing artificial sources
"""

from __future__ import division, print_function

if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')
import h5py
import galsim
import os
import sys
import shutil
import logging
import argparse
import pdb
import photutils

import numpy as np
import matplotlib.pyplot as plt
from pkg_resources import resource_filename
from pickle import dump
from glob import glob
import csv


from astropy.table import Table, Column, vstack
from astropy.io import fits
#from astropy import wcs as astropy_wcs
from fitsio import FITSHDR
import fitsio

from astropy import units
from astropy.coordinates import SkyCoord

from astrometry.libkd.spherematch import match_radec

from tractor.psfex import PsfEx, PsfExModel
from tractor.basics import GaussianMixtureEllipsePSF, RaDecPos

from legacypipe.runbrick import run_brick
from legacypipe.decam import DecamImage
from legacypipe.survey import LegacySurveyData, wcs_for_brick

import obiwan.priors as priors
from obiwan.db_tools import getSrcsInBrick 

from astrometry.util.fits import fits_table, merge_tables
from astrometry.util.ttime import Time

from theValidator.catalogues import CatalogueFuncs

def write_dict(fn,d):
    '''d -- dictionary'''
    w = csv.writer(open(fn, "w"))
    for key, val in d.items():
        w.writerow([key, val])

def read_dict(fn):
    d = {}
    for key, val in csv.reader(open(fn)):
        d[key] = val
    return d

def imshow_stamp(stamp,fn='test.png',galsimobj=True):
    if galsimobj:
        img = stamp.array.copy()
    else:
        img= stamp.copy()
    img=img + abs(img.min())+1 
    plt.imsave(fn,np.log10(img),origin='lower',cmap='gray')
    #plt.imshow(np.log10(img),origin='lower',cmap='gray')
    #plt.savefig(fn)
    #plt.close()
    #print('Wrote %s' % fn)

def plot_radial_profs(fn,profs):
    assert(profs.shape[1] == 3)
    r=np.arange(profs.shape[0])
    for i,lab in zip(range(3),['src','srcnoise','srcnoiseimg']):
        plt.plot(r,profs[:,i],label=lab)
    plt.legend(loc='lower right')
    plt.savefig(fn)
    plt.close()


def ptime(text,t0):
    '''Timer'''    
    tnow=Time()
    print('TIMING:%s ' % text,tnow-t0)
    return tnow



def get_savedir(**kwargs):
  if kwargs['do_skipids'] == 'no':
    return os.path.join(kwargs['decals_sim_dir'],kwargs['objtype'],\
                        kwargs['brickname'][:3], kwargs['brickname'],\
                        "rs%d" % kwargs['rowst'])    
  elif kwargs['do_skipids'] == 'yes':
    return os.path.join(kwargs['decals_sim_dir'],kwargs['objtype'],\
                        kwargs['brickname'][:3], kwargs['brickname'],\
                        "skip_rs%d" % kwargs['rowst'])    

def get_skip_ids(decals_sim_dir, brickname, objtype):
  fns= glob(os.path.join(decals_sim_dir, objtype,
                         brickname[:3], brickname,
                         '*','obiwan','skippedids-*.fits'))
  if len(fns) == 0:
    raise ValueError("no skippedids.fits files exist for this brick %s" % brickname)
  T= CatalogueFuncs().stack(fns, textfile=False)
  return T.ids.astype(str)

def get_fnsuffix(**kwargs):
    return '-{}-{}.fits'.format(kwargs['objtype'], kwargs['brickname'])
                                   #'rs%d' % kwargs['rowst'])

class SimDecals(LegacySurveyData):
    """Top level object that specifying which data to run through pipeline
    
    Same behavior as legacypipe.runs.Dr3DecalsSurvey which chooses which
      CCDs to include. But this also stores all the relevant obiwan
      objects
    
    Args:
      dataset: see definition in 
        https://github.com/legacysurvey/obiwan/blob/master/py/obiwan/test/end_to_end/README.md 
      survey_dir: as used by legacypipe.runbrick.run_brick()
        Defaults to $LEGACY_SURVEY_DIR environment variable.  Where to look for
        files including calibration files, tables of CCDs and bricks, image data
      metacat: fits_table 
        configuration-like params for the simulated sources
      simcat: fits_table
        simulated source catalog for a given brick (not CCD).
      output_dir: legacypipe's outdir
      add_sim_noise: add Poisson noise from the simulated source to the image 
      folding_threshold: how close the simulated source flux is to the requested flux
        make smaller to increase simulated source flux/requested flux
      image_eq_model: referred to as 'testA'
        wherever add a simulated source, replace both image and invvar of the image
        with that of the simulated source only

    Attributes:
      DR: see above 
      metacat: fits_table 
        configuration-like params for the simulated sources
      simcat: fits_table
        simulated source catalog for a given brick (not CCD).
      output_dir: legacypipe's outdir
      add_sim_noise: add Poisson noise from the simulated source to the image 
      folding_threshold: how close the simulated source flux is to the requested flux
        make smaller to increase simulated source flux/requested flux
      image_eq_model: referred to as 'testA'
        wherever add a simulated source, replace both image and invvar of the image
        with that of the simulated source only
    """
    
    def __init__(self, dataset=None, survey_dir=None, metacat=None, simcat=None, output_dir=None,\
		 add_sim_noise=False, folding_threshold=1.e-5, image_eq_model=False):
        super(SimDecals, self).__init__(survey_dir=survey_dir, output_dir=output_dir)
	self.dataset= dataset
	self.metacat = metacat
	self.simcat = simcat
	# Additional options from command line
        self.add_sim_noise= add_sim_noise
        self.folding_threshold= folding_threshold
	self.image_eq_model= image_eq_model
	print('SimDecals: self.image_eq_model=',self.image_eq_model)
        
    def get_image_object(self, t):
	return SimImage(self, t)
    
    #######
    # see legacypipe/runs.py
    def filter_ccds_files(self, fns):
	if self.dataset == 'DR3':
	    return [fn for fn in fns if
		    ('survey-ccds-decals.fits.gz' in fn  or
		     'survey-ccds-nondecals.fits.gz' in fn or
		     'survey-ccds-extra.fits.gz' in fn)]
	#elif self.dataset == 'DR4':
        #	return [fn for fn in fns if
        #			('survey-ccds-dr4-90prime.fits.gz' in fn or
        #			'survey-ccds-dr4-mzlsv2.fits.gz' in fn)]
	elif self.dataset == 'DR5':
	    return fns
    
    def ccds_for_fitting(self, brick, ccds):
	if self.dataset in ['DR3','DR5']:
	    return np.flatnonzero(ccds.camera == 'decam')
	#elif self.dataset == 'DR4':
        #	return np.flatnonzero(np.logical_or(ccds.camera == 'mosaic',
	#						  ccds.camera == '90prime'))
    
    def filter_ccd_kd_files(self, fns):
	"""introduced in DR5"""
	return []
    
def get_srcimg_invvar(stamp_ivar,img_ivar):
    '''stamp_ivar, img_ivar -- galsim Image objects'''
    # Use img_ivar when stamp_ivar == 0, both otherwise
    use_img_ivar= np.ones(img_ivar.array.shape).astype(bool)
    use_img_ivar[ stamp_ivar.array > 0 ] = False
    # First compute using both
    ivar= np.power(stamp_ivar.array.copy(), -1) + np.power(img_ivar.array.copy(), -1) 
    ivar= np.power(ivar,-1) 
    keep= np.ones(ivar.shape).astype(bool)
    keep[ (stamp_ivar.array > 0)*\
          (img_ivar.array > 0) ] = False
    ivar[keep] = 0.
    # Now use img_ivar only where need to
    ivar[ use_img_ivar ] = img_ivar.array.copy()[ use_img_ivar ]
    # return 
    obj_ivar = stamp_ivar.copy()
    obj_ivar.fill(0.)
    obj_ivar+= ivar
    return obj_ivar

def saturation_e(camera):
	# Saturation limit
	d=dict(decam=3e4) # e-
	return d[camera]

def ivar_to_var(ivar,nano2e=None,camera='decam'):
	assert(nano2e is not None)
	flag= ivar == 0.
	var= np.power(ivar, -1)
	# Set 0 ivar pixels to satuation limit
	# var * nano2e^2 = e-^2
	sat= saturation_e(camera) / nano2e**2
	var[flag]= sat
	return var 

class SimImage(DecamImage):
	"""Adds simulated sources to a single exposure

	Similar behavior as legacypipe.decam.DecamImage. Instead of 
		loading images specifically from DECam, this  loads images
		with simulated sources added in 

	Args:
		survey: SimDecals() object
		t: as used by DecamImage
			a single row fits_table for a specific CCD

	Attributes:
		inherits: DecamImage
		t: as used by DecamImage
			a single row fits_table for a specific CCD
	"""

	def __init__(self, survey, t):
		super(SimImage, self).__init__(survey, t)
		self.t = t
                if self.survey.dataset in ['DR3','DR3_eBOSS']:
                    assert('arawgain' in self.t.get_columns())
                    assert(not 'gain' in self.t.get_columns())
                    self.t.rename('arawgain', 'gain')
                elif self.survey.dataset in ['DR5']:
                    assert 'gain' in self.t.get_columns()
                    

	def get_tractor_image(self, **kwargs):
		tim = super(SimImage, self).get_tractor_image(**kwargs)
		if tim is None: # this can be None when the edge of a CCD overlaps
			return tim

		# Seed
		#if 'SEED' in self.survey.metacat.columns:
		#    seed = self.survey.metacat['SEED']
		#else:
		#    seed = None

		objtype = self.survey.metacat.get('objtype')[0]
		objstamp = BuildStamp(tim, gain=self.t.gain, \
				      folding_threshold=self.survey.folding_threshold,\
				      stamp_size= self.survey.metacat.stamp_size)

		# Grab the data and inverse variance images [nanomaggies!]
		tim_image = galsim.Image(tim.getImage())
		tim_invvar = galsim.Image(tim.getInvvar())
		tim_dq = galsim.Image(tim.dq)
		# Also store galaxy sims and sims invvar
		sims_image = tim_image.copy() 
		sims_image.fill(0.0)
		sims_ivar = sims_image.copy()
		# To make cutout for deeplearning
		tim.sims_xy = np.zeros((len(self.survey.simcat),4))-1 
		tim.sims_xyc = np.zeros((len(self.survey.simcat),2))-1
		tim.sims_id = np.zeros(len(self.survey.simcat)).astype(np.int32)-1
		tim.sims_added_flux = np.zeros(len(self.survey.simcat)).astype(float)-1

		# Store simulated galaxy images in tim object 
		# Loop on each object.
		for ii, obj in enumerate(self.survey.simcat):
			# Print timing
			t0= Time()
			if objtype in ['lrg','elg']:
				strin= 'Drawing 1 %s: sersicn=%.2f, rhalf=%.2f, ba=%.2f, phi=%.2f' % \
						(objtype.upper(), obj.sersicn,obj.rhalf,obj.ba,obj.phi)
				print(strin)
			# Before drawing we can check if the obj is near CCD
			#if self.survey.metacat.cutouts[0]: 
			#    #draw_it= isNearCCD(tim,obj,
			#    junk,xx,yy = tim.wcs.wcs.radec2pixelxy(obj.ra,obj.dec)
			#    xx,yy= int(xx),int(yy)
			#    min_stamp_pixels= 16  # 200. / 3600. # arcsec -> deg
			#    obj_bounds= galsim.BoundsI(xmin= xx - min_stamp_pixels/2,\
			#                            xmax= xx + min_stamp_pixels/2,\
			#                            ymin= yy - min_stamp_pixels/2,\
			#                            ymax= yy + min_stamp_pixels/2)
			#    overlap = obj_bounds & tim_image.bounds
			#    # Even the SMALLEST stamp fits entirely within image
			#    # High prob teh full size stamp will  
			#    draw_it= obj_bounds == overlap
			#    #x1, y1 = tim.wcs.positionToPixel(RaDecPos(obj.ra-max_stamp_size/2, obj.dec-max_stamp_size/2))

			if objtype == 'star':
				stamp = objstamp.star(obj)
			elif objtype == 'elg':
				stamp = objstamp.elg(obj)
			elif objtype == 'lrg':
				stamp = objstamp.lrg(obj)
			elif objtype == 'qso':
				stamp = objstamp.qso(obj)
			t0= ptime('Drew the %s' % objtype.upper(),t0)
			#print('I predict we draw it',draw_it)
			# Save radial profiles after draw, addNoise, etc. for unit tests
			#rad_profs=np.zeros((stamp.array.shape[0],3))
			#rad_profs[:,0]= stamp.array.copy()[ stamp.array.shape[0]/2,: ]
			# Want to save flux actually added too
			added_flux= stamp.added_flux
			stamp_nonoise= stamp.copy()
			if self.survey.add_sim_noise:
				#stamp2,stamp3= objstamp.addGaussNoise(stamp, ivarstamp)
				ivarstamp= objstamp.addGaussNoise(stamp)
			# Add source if EVEN 1 pix falls on the CCD
			overlap = stamp.bounds & tim_image.bounds
			add_source = overlap.area() > 0
			# For Deep learning: only add source if entire thing fits on image
			if self.survey.metacat.cutouts[0]:
				# this is a deep learning run
				add_source= stamp.bounds == overlap
			if add_source:
				stamp = stamp[overlap]      
				ivarstamp = ivarstamp[overlap]      
				stamp_nonoise= stamp_nonoise[overlap]
				
				#rad_profs[:,1]= stamp.array.copy()[ stamp.array.shape[0]/2,: ]

				# Zero out invvar where bad pixel mask is flagged (> 0)
				keep = np.ones(tim_dq[overlap].array.shape)
				keep[ tim_dq[overlap].array > 0 ] = 0.
				ivarstamp *= keep
				#tim_invvar[overlap] *= keep # don't modify tim_invvar unless adding stamp ivar

				# Stamp ivar can get messed up at edges
				# especially when needed stamp smaller than args.stamp_size
				cent= int( min(ivarstamp.array.shape)/2 )
    				med= np.median(ivarstamp.array[cent-2:cent+2,cent-2:cent+2].flatten() )
				# 100x median fainter gets majority of star,qso OR elg,lrg profile
				ivarstamp.array[ ivarstamp.array > 100 * med ] = 0.
				# Add stamp to image
				back= tim_image[overlap].copy()
				tim_image[overlap] = back.copy() + stamp.copy()
				# Add variances
				back_ivar= tim_invvar[overlap].copy()
				tot_ivar= get_srcimg_invvar(ivarstamp, back_ivar)
				tim_invvar[overlap] = tot_ivar.copy()

				#rad_profs[:,2]= tim_image[overlap].array.copy()[ stamp.array.shape[0]/2,: ]
				# Save sims info
				tim.sims_xy[ii, :] = [overlap.xmin-1, overlap.xmax-1,
									  overlap.ymin-1, overlap.ymax-1] # galsim 1st index is 1
				tim.sims_xyc[ii, :] = [overlap.trueCenter().x-1, overlap.trueCenter().y-1]
				#tim.sims_radec[ii, :] = [obj.ra,obj.dec]
				tim.sims_id[ii] = obj.id
				tim.sims_added_flux[ii] = added_flux

				# For cutouts we only care about src, background, var (not ivar)
				if self.survey.metacat.cutouts[0]:
					# Data for training: src+noise (cutout) and backgrn (cutout,var,badpix)
					data= np.zeros((stamp.array.shape[0],stamp.array.shape[1],4))
					# FIX ME, add extra rotations for galaxies?
					data[:,:,0]= stamp.array.copy() # src+noise
					#data[:,:,1]= np.sqrt( np.power(stamp.array.copy(),2) )#src+noise var  #ivarstamp.array.copy() 
					data[:,:,1]= back.array.copy() # back
					data[:,:,2]= tim_dq[overlap].array.copy() # bad pix
					data[:,:,3]= stamp_nonoise.array.copy() # Stamp w/out noise, sanity check
					#data[:,:,2]= ivar_to_var(back_ivar.array.copy(),nano2e=objstamp.nano2e) # back var
					#data[:,:,3]= tim_image[overlap].array.copy() # src+noise+background
					#data[:,:,4]= tim_invvar[overlap].array.copy() # src+noise+background_ nvvar
					# Save fn
					brick= os.path.basename(os.path.dirname(self.survey.output_dir))
					hdf5_fn= '%s_%s.hdf5' % (objtype,brick)  #'%s_%d_%s' % (tim.band,obj.id,expid)
					hdf5_fn= os.path.join(self.survey.output_dir,hdf5_fn)
					expid=str(tim.imobj).strip().replace(' ','')
					node= '%s/%s/%s' % (obj.id,tim.band,expid)
					fobj = h5py.File(hdf5_fn, "a")
					dset = fobj.create_dataset(node, data=data,chunks=True)
					for name,val,dtype in zip(\
							['id','flux_added'],\
							[obj.id,added_flux],\
							[np.int32,np.float32]):
						dset.attrs.create(name,val,dtype=dtype)
					#if objtype in ['lrg','elg']:
					#    for name,val in zip(\
					#            ['rhalf','sersicn','phi','ba'],\
					#            [obj.rhalf,obj.sersicn,obj.phi,obj.ba]):
					#        dset.attrs.create(name,val,dtype=np.float32)
						#d.update(dict(rhalf=obj.rhalf,\
						#              sersicn=obj.sersicn,\
						#              phi=obj.phi,\
						#              ba=obj.ba))
					print('Saved %s to %s' % (node,hdf5_fn))
					#np.save(fn+'.npy',data,allow_pickle=False)
					# Save enough metadata to classify image quality later
					#x1,x2,y1,y2= tuple(tim.sims_xy[ii,:])
					#xc,yc= tuple(tim.sims_xyc[ii,:])
					#d = dict(band=tim.band,\
					#		 expid=expid,\
					#         addedflux= added_flux,\
					#         id=obj.id,\
					#         ra=obj.ra,\
					#         dec=obj.dec)
							 #xc=xc,yc=yc,\
							 #(x1=x1,x2=x2,y1=y1,y2=y2,\
							 #gflux=obj.gflux,\
							 #rflux=obj.rflux,\
							 #zflux=obj.zflux)
					#write_dict(fn+'.csv',d)
					# Write sanity checks if they don't exists
					#fns= glob(os.path.join(self.survey.output_dir,'*_src.fits'))
					#if len(fns) == 0:
					#    # Also write fits file for easier image stretching
					#    fitsio.write(fn+'_src.fits',data[...,0],clobber=True)
					#    fitsio.write(fn+'_src_invvar.fits',data[...,1],clobber=True)
					#    fitsio.write(fn+'_img.fits',data[...,2],clobber=True)
					#    fitsio.write(fn+'_img_invvar.fits',data[...,3],clobber=True)
					#    fitsio.write(fn+'_srcimg.fits',data[...,4],clobber=True)
					#    fitsio.write(fn+'_srcimg_invvar.fits',data[...,5],clobber=True)
					#    # Draw Radial Profiles
					#    plot_radial_profs(fn+'_profiles.png',rad_profs)
				
				#Extra
				sims_image[overlap] += stamp.copy() 
				sims_ivar[overlap] += ivarstamp.copy()
				
					
				#print('HACK!!!')
				#galsim.fits.write(stamp, 'stamp-{:02d}.fits'.format(ii), clobber=True)
				#galsim.fits.write(ivarstamp, 'ivarstamp-{:02d}.fits'.format(ii), clobber=True)

				if np.min(sims_ivar.array) < 0:
					log.warning('Negative invvar!')
					import pdb ; pdb.set_trace()
		tim.sims_image = sims_image.array
		tim.sims_inverr = np.sqrt(sims_ivar.array)
		tim.sims_xy = tim.sims_xy.astype(int)
		tim.sims_xyc = tim.sims_xyc.astype(int)
		# Can set image=model, ivar=1/model for testing
		if self.survey.image_eq_model:
			tim.data = sims_image.array.copy()
			tim.inverr = np.zeros(tim.data.shape)
			tim.inverr[sims_image.array > 0.] = np.sqrt(1./sims_image.array.copy()[sims_image.array > 0.]) 
		else:
			tim.data = tim_image.array
			tim.inverr = np.sqrt(tim_invvar.array)
		 
		#print('HACK!!!')
		#galsim.fits.write(invvar, 'invvar.fits'.format(ii), clobber=True)
		#import pdb ; pdb.set_trace()
		return tim

class BuildStamp():
	"""Does the drawing of simulated sources on a single exposure

	Args: 
		tim: Tractor Image Object for a specific CCD
		gain: gain of the CCD
		folding_threshold: how close the simulated source flux is to the requested flux
			make smaller to increase simulated source flux/requested flux
		stamp_size: pixels, width and height of simulated images

	Attributes:
		band: g,r,z
		stamp_size: pixels, width and height of simulated images
		gsparams: galsim object that configures how accurate simulated source will be
		gsdeviate: galsim object that configures its random number generator
		wcs: WCS from tim
		psf: psf from tim
		galsim_wcs: wcs repackaged into galsim compatible object
		zpscale: conversion factor 'nanomaggies' to 'ADU'
		nano2e: conversion factor 'nanomaggies' to 'e-'
	"""

	def __init__(self,tim, gain=4.0, folding_threshold=1.e-5, stamp_size=None):
		self.band = tim.band.strip()
		self.stamp_size = stamp_size
		# GSParams should be used when galsim object is initialized
		# MAX size for sersic: 
		# https://github.com/GalSim-developers/GalSim/pull/450/commits/755bcfdca25afe42cccfd6a7f8660da5ecda2a65
		MAX_FFT_SIZE=1048576L #2^16=65536
		self.gsparams = galsim.GSParams(maximum_fft_size=MAX_FFT_SIZE,\
										folding_threshold=folding_threshold) 
		#print('FIX ME!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
		self.gsdeviate = galsim.BaseDeviate()
		#if seed is None:
		#    self.gsdeviate = galsim.BaseDeviate()
		#else:
		#    self.gsdeviate = galsim.BaseDeviate(seed)

                self.wcs = tim.getWcs()
		self.psf = tim.getPsf()
		# Tractor wcs object -> galsim wcs object
		temp_hdr = FITSHDR()
		subwcs = tim.wcs.wcs.get_subimage(tim.wcs.x0, tim.wcs.y0,
								  int(tim.wcs.wcs.get_width())-tim.wcs.x0,
								  int(tim.wcs.wcs.get_height())-tim.wcs.y0)
		subwcs.add_to_header(temp_hdr)
		# Galsim uses astropy header, not fitsio
		hdr = fits.Header()
		for key in temp_hdr.keys(): hdr[key]=temp_hdr[key]
		self.galsim_wcs = galsim.GSFitsWCS(header=hdr)
		del subwcs,temp_hdr,hdr
		
		# zpscale equivalent to magzpt = self.t.ccdzpt+2.5*np.log10(self.t.exptime)
		self.zpscale = tim.zpscale      # nanomaggies-->ADU conversion factor
		self.nano2e = self.zpscale*gain # nanomaggies-->electrons conversion factor

	def setlocal(self,obj):
		"""Get the pixel positions, local wcs, local PSF.""" 

		xx, yy = self.wcs.positionToPixel(RaDecPos(obj.get('ra'), obj.get('dec')))
		self.pos = galsim.PositionD(xx, yy)
		self.xpos = int(self.pos.x)
		self.ypos = int(self.pos.y)
		self.offset = galsim.PositionD(self.pos.x-self.xpos, self.pos.y-self.ypos)

		# galsim.drawImage() requires local (linear) wcs
		self.localwcs = self.galsim_wcs.local(image_pos=self.pos)
		#cd = self.wcs.cdAtPixel(self.pos.x, self.pos.y)
		#self.pixscale = np.sqrt(np.linalg.det(cd))*3600.0
		
		# Get the local PSF
		psfim = self.psf.getPointSourcePatch(self.xpos, self.ypos).getImage()
		#plt.imshow(psfim) ; plt.show()
		
		# make galsim PSF object
		self.localpsf = galsim.InterpolatedImage(galsim.Image(psfim), wcs=self.galsim_wcs,\
												 gsparams=self.gsparams)

	def addGaussNoise(self, stamp):
		"""Adds gaussian noise to perfect source (in place)

                STAMP and IVARSTAMP are in units of nanomaggies and 
                  1/nanomaggies**2, respectively.
		
                Returns:
                  invvar for the stamp
		"""
		#stamp= stamp_backup.copy()
		#ivarstamp= ivarstamp_backup.copy()

		
		#varstamp = ivarstamp.copy()
		#ivarstamp.invertSelf() # input data, convert to variance
		#ivarstamp *= self.nano2e**2 # [electron^2]
			 
		# Add the variance of the object to the variance image (in electrons).
		stamp *= self.nano2e       # [noiseless stamp, electron]
		stamp_var = galsim.Image(np.sqrt(stamp.array**2), wcs=self.galsim_wcs) 
		stamp_var.setOrigin(galsim.PositionI(stamp.xmin, stamp.ymin))

		# Add Poisson noise
		noise = galsim.VariableGaussianNoise(self.gsdeviate, stamp_var)
		#stamp2= stamp.copy()
		stamp.addNoise(noise)
		#stamp3= stamp2.copy()
		#c=np.random.normal(loc=0,scale=np.sqrt(objvar.array),size=objvar.array.shape)
		#noise = galsim.Image(c, wcs=self.galsim_wcs)
		#noise.setOrigin(galsim.PositionI(stamp.xmin, stamp.ymin))
		#stamp3+= noise
		
		# Variance of stamp+noise
		stamp_var = stamp.copy()
		stamp_var.fill(0.)
		stamp_var+= np.abs( stamp.array.copy() )
		
		#imshow_stamp(stamp,fn='std.png')
		#imshow_stamp(stamp_backup,'img.png')
		
		#b = galsim.Image(np.zeros(stamp.array.shape), wcs=self.galsim_wcs) 
		#b.array+= stamp.array.copy() 
		#b.array+= stamp_backup.array.copy() 
		#b= stamp.array.copy() + stamp_backup.array.copy()
		#b= stamp.copy()
		#b.drawImage(stamp_backup.copy(),add_to_image=True)
		#imshow_stamp(b,fn='std_img.png')
		# hists
		#for data,nam in zip([stamp.array.copy(),stamp_backup.array.copy(),b],['std','img','std_img']):
		#    j=plt.hist(data)
		#    plt.savefig(nam+'_hist.png')
		#    plt.close(nam+'_hist.png')
		# Convert back to [nanomaggies]
		stamp /= self.nano2e      
		#stamp2 /= self.nano2e      
		#stamp3 /= self.nano2e      
		stamp_var /= self.nano2e**2

		#ivarstamp = varstamp.copy()
		stamp_var.invertSelf()
		# Remask pixels that were masked in the original inverse variance stamp.
		#ivarstamp *= mask
		# This is now inv variance
		return stamp_var

	def convolve_and_draw(self,obj):
		"""Convolve the object with the PSF and then draw it."""
		obj = galsim.Convolve([obj, self.localpsf], gsparams=self.gsparams)
		# drawImage() requires local wcs
		#try:
		if self.stamp_size is None:
			stamp = obj.drawImage(offset=self.offset, wcs=self.localwcs,method='no_pixel')
		else:
			stamp = obj.drawImage(offset=self.offset, wcs=self.localwcs,method='no_pixel',\
								  nx=self.stamp_size,ny=self.stamp_size)
		
		#except SystemExit:
		#except BaseException:
		#    #logging.error(traceback.format_exc())
		#    print('got back drawImage!')
		#    raise ValueError
		#try: 
		#except:
		#    print("Unexpected error:", sys.exc_info()[0])
		#    raise
		stamp.setCenter(self.xpos, self.ypos)
		return stamp

	def star(self,obj):
		"""Render a star (PSF)."""
		log = logging.getLogger('decals_sim')
		# Use input flux as the 7'' aperture flux
		self.setlocal(obj)
		psf = self.localpsf.withFlux(1.)
		if self.stamp_size is None:
			stamp = psf.drawImage(offset=self.offset, wcs=self.localwcs, method='no_pixel')
		else:
			stamp = psf.drawImage(offset=self.offset, wcs=self.localwcs, method='no_pixel',\
								  nx=self.stamp_size,ny=self.stamp_size)
		# Fraction flux in 7'', FIXED pixelscale
		diam = 7/0.262
		# Aperture fits on stamp
		width= stamp.bounds.xmax-stamp.bounds.xmin
		height= stamp.bounds.ymax-stamp.bounds.ymin
		if diam > width and diam > height:
			nxy= int(diam)+2
			stamp = psf.drawImage(nx=nxy,ny=nxy, offset=self.offset, wcs=self.localwcs, method='no_pixel')
		assert(diam <= float(stamp.bounds.xmax-stamp.bounds.xmin))
		assert(diam <= float(stamp.bounds.ymax-stamp.bounds.ymin))
		# Aperture photometry
		apers= photutils.CircularAperture((stamp.trueCenter().x,stamp.trueCenter().y), r=diam/2)
		apy_table = photutils.aperture_photometry(stamp.array, apers)
		apflux= np.array(apy_table['aperture_sum'])[0]
		# Incrase flux so input flux contained in aperture
		flux = obj.get(self.band+'flux')*(2.-apflux/stamp.added_flux) # [nanomaggies]
		psf = self.localpsf.withFlux(flux)
		if self.stamp_size is None:
			stamp = psf.drawImage(offset=self.offset, wcs=self.localwcs, method='no_pixel')
		else:
			stamp = psf.drawImage(offset=self.offset, wcs=self.localwcs, method='no_pixel',\
								  nx=self.stamp_size,ny=self.stamp_size)
		# stamp looses less than 0.01% of requested flux
		if stamp.added_flux/flux <= 0.9999:
			log.warning('stamp lost more than 0.01 percent of requested flux, stamp_flux/flux=%.7f',stamp.added_flux/flux)
		# test if obj[self.band+'FLUX'] really is in the 7'' aperture
		#apers= photutils.CircularAperture((stamp.trueCenter().x,stamp.trueCenter().y), r=diam/2)
		#apy_table = photutils.aperture_photometry(stamp.array, apers)
		#apflux= np.array(apy_table['aperture_sum'])[0]
		#print("7'' flux/input flux= ",apflux/obj[self.band+'FLUX'])
		
		# Convert stamp's center to its corresponding center on full tractor image
		stamp.setCenter(self.xpos, self.ypos)
		return stamp 

	def elg(self,obj):
		"""Create an ELG (disk-like) galaxy."""
		# Create localpsf object
		self.setlocal(obj)
		objflux = obj.get(self.band+'flux') # [nanomaggies]
		try:
			galobj = galsim.Sersic(float(obj.get('sersicn')), half_light_radius=float(obj.get('rhalf')),\
								flux=objflux, gsparams=self.gsparams)
		except:
			raise ValueError 
		galobj = galobj.shear(q=float(obj.get('ba')), beta=float(obj.get('phi'))*galsim.degrees)
		stamp = self.convolve_and_draw(galobj)
		return stamp

	def lrg(self,obj):
		"""Create an LRG just like did for ELG"""
		return self.elg(obj)

	def qso(self,obj):
		"""Create a QSO just like a star"""
		return self.star(obj)



def flag_nearest_neighbors(Samp, radius_in_deg=5./3600):
  """Returns Sample indices to keep (have > dist separations) and indices to skip

  Returns:
    tuple: keep,skip: indices of Samp to keep and skip
  """
  flag_set=set()
  all_indices= range(len(Samp))
  for cnt in all_indices:
      if cnt in flag_set:
          continue
      else:
          I,J,d = match_radec(Samp.ra[cnt],Samp.dec[cnt],
                              Samp.ra,Samp.dec, 5./3600,
                              notself=False,nearest=False)
          # Remove all Samp matches (J), minus the reference ra,dec
          flag_inds= set(J).difference(set( [cnt] ))
          if len(flag_inds) > 0:
              flag_set= flag_set.union(flag_inds)
  keep= list( set(all_indices).difference(flag_set) )
  return keep, list(flag_set)

#def build_simcat(nobj=None, brickname=None, brickwcs=None, meta=None, seed=None, noOverlap=True):
def build_simcat(Samp=None,brickwcs=None, meta=None):
    """Creates the simulated source catalog for a given brick (not CCD).

	The WCS for the brick (not CCD) is used to convert ra,dec of source
		to x,y pixel location in brickspace

	Args:
		Samp: fits_table for the properties of sources in the brick
			usually a subset of all sources in the brick determined by
			rowstart (rs)
		brickwcs: WCS object for the brick
		meta: 'metacat' table 
			fits_table with configuration-like params for the simulated sources

	Returns: 
		tuple of
		cat:
		skipping_ids:
	"""
    log = logging.getLogger('decals_sim')

    #rand = np.random.RandomState(seed)

    # Assign central coordinates uniformly but remove simulated sources which
    # are too near to one another.  Iterate until we have the requisite number
    # of objects.
    #bounds = brickwcs.radec_bounds()
    #ra = rand.uniform(bounds[0], bounds[1], nobj)
    #dec = rand.uniform(bounds[2], bounds[3], nobj)
    i_keep,i_skip= flag_nearest_neighbors(Samp, radius_in_deg=5./3600)
    skipping_ids= Samp.get('id')[i_skip]
    log.info('sources %d, keeping %d, flagged as nearby %d' % (len(Samp),len(i_keep),len(i_skip)))
    Samp.cut(i_keep)
    
    xxyy = brickwcs.radec2pixelxy(Samp.ra,Samp.dec)

    #cat = Table()
    #cat['ID'] = Column(Samp.get('id'),dtype='i4') #np.arange(nobj, dtype='i4'))
    #cat['RA'] = Column(Samp.ra, dtype='f8')
    #cat['DEC'] = Column(Samp.dec, dtype='f8')
    #cat['X'] = Column(xxyy[1][:], dtype='f4')
    #cat['Y'] = Column(xxyy[2][:], dtype='f4')
    cat = fits_table()
    for key in ['id','seed','ra','dec']:
        cat.set(key, Samp.get(key))
    cat.set('x', xxyy[1][:])
    cat.set('y', xxyy[2][:])

    typ=meta.get('objtype')[0]
    # Mags
    for key in ['g','r','z']:
        if (meta.bright_galaxies[0]) and (typ in ['elg','lrg']):
            # Galaxies get stellar flux, so they are easy to see in images
            cat.set('%sflux' % key, 1E9*10**(-0.4*Samp.get('%s_%s' % ('star',key))) ) # [nanomaggies]
        else:
            # Actual galaxy colors
            cat.set('%sflux' % key, 1E9*10**(-0.4*Samp.get('%s_%s' % (typ,key))) ) # [nanomaggies]
        # Galaxy Properties
    if typ in ['elg','lrg']:
        for key,tab_key in zip(['sersicn','rhalf','ba','phi'],['n','re','ba','pa']):
            cat.set(key, Samp.get('%s_%s'%(typ,tab_key) ))
        # Sersic n: GALSIM n = [0.3,6.2] for numerical stability,see
        # https://github.com/GalSim-developers/GalSim/issues/{325,450}
        # I'll use [0.4,6.1]
        vals= cat.sersicn
        vals[cat.sersicn < 0.4] = 0.4
        vals[cat.sersicn > 6.1] = 6.1
        cat.set('sersicn',vals)
        #cat['R50_1'] = Column(Samp.rhalf, dtype='f4')
        #cat['BA_1'] = Column(Samp.ba, dtype='f4')
        #cat['PHI_1'] = Column(Samp.phi, dtype='f4')
    return cat, skipping_ids



def get_parser():
    '''return parser object, tells it what options to look for
    options can come from a list of strings or command line'''
    parser = argparse.ArgumentParser(formatter_class=argparse.
                                     ArgumentDefaultsHelpFormatter,
                                     description='DECaLS simulations.')
    parser.add_argument('--dataset', type=str, choices=['DR5','DR3', 'DR3_eBOSS'], required=True, help='see definitions in obiwan/test/README.md') 
    parser.add_argument('-o', '--objtype', type=str, choices=['star','elg', 'lrg', 'qso'], default='star', required=True) 
    parser.add_argument('-b', '--brick', type=str, default='2428p117', required=True)
    parser.add_argument('--outdir', default='./', required=False)
    parser.add_argument('-n', '--nobj', type=int, default=500, metavar='', 
                        help='number of objects to simulate (required input)')
    parser.add_argument('-rs', '--rowstart', type=int, default=0, metavar='', 
                        help='zero indexed, row of ra,dec,mags table, after it is cut to brick, to start on')
    parser.add_argument('--do_skipids', type=str, choices=['no','yes'],default='no', help='inject skipped ids for brick, otherwise run as usual')
    parser.add_argument('--randoms_db', default='obiwan_elg', help='desi db table name for randoms')
    parser.add_argument('--randoms_from_fits', default=None, help='set to read randoms from fits file instead of scidb2.nersc.gov db, set to absolute path of local fits file on computer')
    parser.add_argument('--prefix', type=str, default='', metavar='', 
                        help='tells which input sample to use')
    #parser.add_argument('-ic', '--ith_chunk', type=long, default=None, metavar='', 
    #                    help='run the ith chunk, 0-999')
    #parser.add_argument('-c', '--nchunk', type=long, default=1, metavar='', 
    #                    help='run chunks 0 to nchunk')
    parser.add_argument('-t', '--threads', type=int, default=1, metavar='', 
                        help='number of threads to use when calling The Tractor')
    #parser.add_argument('-s', '--seed', type=long, default=None, metavar='', 
    #                    help='random number seed, determines chunk seeds 0-999')
    parser.add_argument('-z', '--zoom', nargs=4, default=(0, 3600, 0, 3600), type=int, metavar='', 
                        help='see runbrick.py; (default is 0 3600 0 3600)')
    parser.add_argument('-survey-dir', '--survey_dir', metavar='', 
                        help='Location of survey-ccds*.fits.gz')
    #parser.add_argument('--rmag-range', nargs=2, type=float, default=(18, 26), metavar='', 
    #                    help='r-band magnitude range')
    parser.add_argument('--add_sim_noise', action="store_true", help="set to add noise to simulated sources")
    parser.add_argument('--folding_threshold', type=float,default=1.e-5,action="store", help="for galsim.GSParams")
    parser.add_argument('-testA','--image_eq_model', action="store_true", help="set to set image,inverr by model only (ignore real image,invvar)")
    parser.add_argument('--all-blobs', action='store_true', 
                        help='Process all the blobs, not just those that contain simulated sources.')
    parser.add_argument('--stage', choices=['tims', 'image_coadds', 'srcs', 'fitblobs', 'coadds'],
                        type=str, default=None, metavar='', help='Run up to the given stage')
    parser.add_argument('--early_coadds', action='store_true',default=False,
                        help='add this option to make the JPGs before detection/model fitting')
    parser.add_argument('--cutouts', action='store_true',default=False,
                        help='Stop after stage tims and save .npy cutouts of every simulated source')
    parser.add_argument('--stamp_size', type=int,action='store',default=64,\
                        help='Stamp/Cutout size in pixels')
    parser.add_argument('--bright_galaxies', action='store_true',
                        help='Galaxies get stellar colors for DEEP Obiwan training')
    parser.add_argument('--bricklist',action='store',default='bricks-eboss-ngc.txt',\
                        help='if using mpi4py, $LEGACY_SURVEY_DIR/bricklist')
    parser.add_argument('--nproc', type=int,action='store',default=1,\
                        help='if using mpi4py')
    parser.add_argument('-v', '--verbose', action='store_true', help='toggle on verbose output')
    return parser
 
def create_metadata(kwargs=None):
    """fits_table with configuration-like params for the simulated sources
    
    TODO: Should metacat table have a rowstart column?
    TODO: One metacat table per brick, instead of one per `rs*` directory?
    
    Args:
      kwargs: configuration-like params for the simulated sources
        {'brickname': which chunk of sky
        'objtype': star,elg,lrg,qso
        'nobj': number of simulated sources for this run
        'stamp_size': pixels, width and height of simulated images
        'cutouts': whether .npy cutouts of every simulated source were written
        'bright_galaxies': whether bright_galaxies flag is set
        }
    
    Returns:
      Nothing
      writes the 'metacat' fits_table to disk and stores it
      in the kwargs input arg
    """
    assert(kwargs is not None)
    log = logging.getLogger('decals_sim')
    # Pack the input parameters into a meta-data table and write out.
    #metacols = [
    #    ('BRICKNAME', 'S10'),
    #    ('OBJTYPE', 'S10'),
    #    ('NOBJ', 'i4'),
    #    ('CHUNKSIZE', 'i2'),
    #    ('NCHUNK', 'i2'),
    #    ('ZOOM', 'i4', (4,)),
    #    ('SEED', 'S20'),
    #    ('RMAG_RANGE', 'f4', (2,))]
    #metacat = Table(np.zeros(1, dtype=metacols))
    metacat = fits_table()
    for key in ['brickname','objtype']: #,'nchunk']:
	metacat.set(key, np.array( [kwargs[key]] ))
	metacat.set('nobj', np.array( [kwargs['args'].nobj] ))
	metacat.set('zoom', np.array( [kwargs['args'].zoom] ))
	metacat.set('cutouts', np.array( [kwargs['args'].cutouts] ))
	metacat.set('stamp_size', np.array( [kwargs['args'].stamp_size] ))
	metacat.set('bright_galaxies', np.array( [kwargs['args'].bright_galaxies] ))
	#metacat['RMAG_RANGE'] = kwargs['args'].rmag_range
	#if not kwargs['args'].seed:
	#    log.info('Random seed = {}'.format(kwargs['args'].seed))
	#    metacat['SEED'] = kwargs['args'].seed
    #metacat_dir = os.path.join(kwargs['decals_sim_dir'], kwargs['objtype'],kwargs['brickname'][:3],kwargs['brickname'])    
    metacat_dir = get_savedir(**kwargs)
    if not os.path.exists(metacat_dir): 
	os.makedirs(metacat_dir)
    
    metafile = os.path.join(metacat_dir, 'metacat'+get_fnsuffix(**kwargs))
    log.info('Writing {}'.format(metafile))
    if os.path.isfile(metafile):
	os.remove(metafile)
    metacat.writeto(metafile)
    
    # Store new stuff
    kwargs['metacat']=metacat
    kwargs['metacat_dir']=metacat_dir


def create_ith_simcat(d=None):
	"""Write 'simcat' and 'skipped_ids' tables for a given sample of sources

        Args:
          d: {'Samp': fits_table for the properties of sources in the brick
	     'brickwcs': WCS object for the brick
	     'metacat': fits_table with configuration params for the simulated sources
             }
        
        Returns:
	  Nothing, saves the 'simcat' and 'skipped_ids' tables
	  Adds 'simcat' table to dict 'd'
	"""
	assert(d is not None)
	log = logging.getLogger('decals_sim')
	#chunksuffix = '{:02d}'.format(ith_chunk)
	# Build and write out the simulated object catalog.
	#seed= d['seeds'][ith_chunk]
	#simcat = build_simcat(d['nobj'], d['brickname'], d['brickwcs'], d['metacat'], seed)
	simcat, skipped_ids = build_simcat(Samp=d['Samp'],brickwcs=d['brickwcs'],meta=d['metacat'])
	# Simcat 
	simcat_dir = get_savedir(**d) #os.path.join(d['metacat_dir'],'row%d-%d' % (rowstart,rowend)) #'%3.3d' % ith_chunk)    
	if not os.path.exists(simcat_dir): 
		os.makedirs(simcat_dir)
	#simcatfile = os.path.join(simcat_dir, 'simcat-{}-{}-row{}-{}.fits'.format(d['brickname'], d['objtype'],rowstart,rowend)) # chunksuffix))
	simcatfile = os.path.join(simcat_dir, 'simcat'+get_fnsuffix(**d))
	if os.path.isfile(simcatfile):
		os.remove(simcatfile)
	simcat.writeto(simcatfile)
	log.info('Wrote {}'.format(simcatfile))
	# Skipped Ids
	if len(skipped_ids) > 0:
		skip_table= fits_table()
		skip_table.set('ids',skipped_ids)
		name= os.path.join(simcat_dir,'skippedids'+get_fnsuffix(**d))
		if os.path.exists(name):
			os.remove(name)
			log.info('Removed %s' % name)
		skip_table.writeto(name)
		log.info('Wrote {}'.format(name))
	# add to dict
	d['simcat']= simcat
	d['simcat_dir']= simcat_dir

def get_runbrick_setup(**kwargs):
    """Convert runbrick.py cmd line options into `**kwargs` for run_brick()
    
    The command line options depend on the Data Release (e.g. the
      legacypipe code version. The cmd line options associated with 
      each DR get modified and repackaged into a dict in 
      legacypipe.runbrick so this converter is required to call run_brick
      appropriately
    
    Args:
      **kwargs: dict of the cmd line options to obiwan.kenobi.py
    
    Returns:
      dict to use when calling legacypipe.runbrick.run_brick like
        run_brick(brickname, survey, `**dict`)		
    """
    dataset= kwargs['dataset']
    assert(dataset in ['DR5','DR3','DR3_eBOSS'])
    from legacypipe.runbrick import get_runbrick_kwargs
    from legacypipe.runbrick import get_parser as get_runbrick_parser
    zm= kwargs['zoom']
    cmd_line= ['--no-write', '--skip','--force-all',
	       '--zoom','%d' % zm[0],'%d' % zm[1],'%d' % zm[2],'%d' % zm[3],
	       '--no-wise', '--threads','%d' % kwargs['threads']]
    if kwargs['early_coadds']:
        cmd_line += ['--early-coadds']
    if kwargs['stage']:
        cmd_line += ['--stage', '%s' % kwargs['stage']]
    if dataset == 'DR3':
        #cmd_line += ['--run', 'dr3', '--hybrid-psf','--nsigma', '6']
        cmd_line += ['--run', 'dr3','--nsigma', '6']
    elif dataset == 'DR5':
        # defaults: rex (use --simp), nsigma 6 ,hybrid-psf (--no-hybrid-psf otherwise)
        # depth cut already done (use --depth-cut to do depth cut anyway)
        cmd_line += ['--run', 'dr5'] 
    
    rb_parser= get_runbrick_parser()
    rb_opt = rb_parser.parse_args(args=cmd_line)
    rb_optdict = vars(rb_opt)
    # remove keys as Dustin' does
    _= rb_optdict.pop('ps', None)
    _= rb_optdict.pop('verbose',None)
    _, rb_kwargs= get_runbrick_kwargs(**rb_optdict)
    return rb_kwargs

def do_one_chunk(d=None):
    """Runs the legacypipe/Tractor pipeline on images with simulated sources
    
    Args:
      d: {'args': obiwan.kenobi.py cmd line argparse.Namespace object
          'brickname': chunk of sky
          'metacat': fits_table	configuration params for the simulated sources
          'simcat': fits_table simulated source catalog for a given brick (not CCD).
    
    Note:
      runb_brick() is 'main' for the legacypipe/Tractor pipeline
    
    Returns:
      Nothing, but this func end ups writing out all the obiwan results 
    """
    assert(d is not None)
    simdecals = SimDecals(dataset=d['args'].dataset,\
			  metacat=d['metacat'], simcat=d['simcat'], output_dir=d['simcat_dir'], \
			  add_sim_noise=d['args'].add_sim_noise, folding_threshold=d['args'].folding_threshold,\
			  image_eq_model=d['args'].image_eq_model)
    # Use Tractor to just process the blobs containing the simulated sources.
    if d['args'].all_blobs:
	blobxy = None
    else:
	blobxy = zip(d['simcat'].get('x'), d['simcat'].get('y'))
    
    # Default runbrick call sequence
    obiwan_kwargs= vars(d['args']) 
    runbrick_kwargs= get_runbrick_setup(**obiwan_kwargs)
    # Obiwan modifications
    runbrick_kwargs.update(blobxy=blobxy)
    #plotbase='obiwan')
    print('Calling run_brick with: ')
    print('brickname= %s' % d['brickname'])
    print('simdecals= ',simdecals)
    print('runbrick_kwards= ',runbrick_kwargs)
    # Run it: run_brick(brick, survey obj, **kwargs)
    run_brick(d['brickname'], simdecals, **runbrick_kwargs)

def dobash(cmd):
    print('UNIX cmd: %s' % cmd)
    if os.system(cmd): raise ValueError

def do_ith_cleanup(d=None):
    """Moves all obiwan+legacypipe outputs to a new directory stucture

    Uses rsync to move everthing, if all the rsync's succeed then all the 
        original files and directories are removed 

    Args:
        d: dict with keys brickname, simcat_dir
    """
    assert(d is not None) 
    log = logging.getLogger('decals_sim')
    log.info('Cleaning up...')
    brick= d['brickname']
    bri= brick[:3]
    outdir= d['simcat_dir']
    # coadd/123/1238p245/* -> coadd/
    # metrics/123/* -> metrics/
    # tractor/123/* -> tractor/
    # tractor-i/123/* -> tractor-i/
    # *.fits -> obiwan/
    try:
        dobash('rsync -av %s/coadd/%s/%s/* %s/coadd/' % 
                (outdir,bri,brick, outdir))
        for name in ['metrics','tractor','tractor-i']:
            dobash('rsync -av %s/%s/%s/* %s/%s/' % 
                    (outdir,name,bri, outdir,name))
        dobash('mkdir -p %s/obiwan' % outdir)
        dobash('rsync -av %s/*.fits %s/obiwan/' % 
                (outdir,outdir))
    except:
        raise ValueError('issue repackaging %s' % output_dir)

    # If here, then safe to remove originals
    dobash('rm -r %s/coadd/%s' % (outdir,bri))
    for name in ['metrics','tractor','tractor-i']:
        dobash('rm -r %s/%s/%s' % (outdir,name,bri))
    dobash('rm %s/*.fits' % outdir)
    # Only "rowstars 0" keeps its coadds
    if d['rowst'] != 0:
        dobash('rm %s/coadd/*.fits.fz' % outdir)
        dobash('rm %s/coadd/*.fits' % outdir)
        


def get_sample_fn(brick,decals_sim_dir,prefix=''):
    fn= os.path.join(decals_sim_dir,'input_sample','bybrick','%ssample_%s.fits' % (prefix,brick))
    return fn
    #return os.path.join(decals_sim_dir,'softlinked_table') #'sample-merged.fits')

def main(args=None):
    """Main routine which parses the optional inputs."""
    t0= Time()
    # Command line options
    if args is None:
        # Read from cmd line
        parser= get_parser()  
        args = parser.parse_args(args=args)
    else:
        # args is already a argparse.Namespace obj
        pass 
    # Print calling sequence
    print('Args:', args)   
    if args.cutouts:
        args.stage = 'tims'
    # Setup loggers
    if args.verbose:
        lvl = logging.DEBUG
    else:
        lvl = logging.INFO
    logging.basicConfig(level=lvl, stream=sys.stdout) #,format='%(message)s')
    log = logging.getLogger('decals_sim')
    # Sort through args 
    log.info('decals_sim.py args={}'.format(args))
    #max_nobj=500
    #max_nchunk=1000
    #if args.ith_chunk is not None: assert(args.ith_chunk <= max_nchunk-1)
    #assert(args.nchunk <= max_nchunk)
    #assert(args.nobj <= max_nobj)
    #if args.ith_chunk is not None: 
    #    assert(args.nchunk == 1) #if choose a chunk, only doing 1 chunk
    if args.nobj is None:
        parser.print_help()
        sys.exit(1)

    brickname = args.brick
    objtype = args.objtype
    maxobjs = args.nobj

    for obj in ('LSB'):
        if objtype == obj:
            log.warning('{} objtype not yet supported!'.format(objtype))
            return 0

    # Output dir
    decals_sim_dir = args.outdir
        
    #nchunk = args.nchunk
    #rand = np.random.RandomState(args.seed) # determines seed for all chunks
    #seeds = rand.random_integers(0,2**18, max_nchunk)

    log.info('Object type = {}'.format(objtype))
    #log.info('Number of objects = {}'.format(nobj))
    #log.info('Number of chunks = {}'.format(nchunk))

    # Optionally zoom into a portion of the brick
    survey = LegacySurveyData()
    brickinfo = survey.get_brick_by_name(brickname)
    print(brickname)
    brickwcs = wcs_for_brick(brickinfo)
    W, H, pixscale = brickwcs.get_width(), brickwcs.get_height(), brickwcs.pixel_scale()

    log.info('Brick = {}'.format(brickname))
    if args.zoom is not None: # See also runbrick.stage_tims()
        (x0, x1, y0, y1) = args.zoom
        W = x1 - x0
        H = y1 - y0
        brickwcs = brickwcs.get_subimage(x0, y0, W, H)
        log.info('Zoom (pixel boundaries) = {}'.format(args.zoom))
    targetrd = np.array([brickwcs.pixelxy2radec(x, y) for x, y in
                         [(1,1), (W,1), (W,H), (1,H), (1,1)]])

    radec_center = brickwcs.radec_center()
    log.info('RA, Dec center = {}'.format(radec_center))
    log.info('Brick = {}'.format(brickname))
    t0= ptime('First part of Main()',t0)
    
    #if args.ith_chunk is not None: 
    #    chunk_list= [args.ith_chunk]
    #else: 
    #    chunk_list= range(nchunk)
    #chunk_list= [ int((args.rowstart)/maxobjs) ]

    # Ra,dec,mag table
    print('before PSQL')
    if args.randoms_from_fits:
        # Non PSQL way
        #fn= get_sample_fn(brickname,decals_sim_dir,prefix=args.prefix)
        #fn=os.path.join(decals_sim_dir,
        #                'elg_randoms/rank_1_seed_1.fits') 
        Samp= fits_table(args.randoms_from_fits)
        Samp.rename('%s_rhalf' % objtype,'%s_re' % objtype)
    else:
      if args.do_skipids == 'no':
        Samp= getSrcsInBrick(brickname,objtype, db_table=args.randoms_db)
      elif args.do_skipids == 'yes':
        skip_ids= get_skip_ids(decals_sim_dir, brickname, objtype)
        Samp= getSrcsInBrick(brickname,objtype, db_table=args.randoms_db,
                             skipped_ids= skip_ids)
    print('after PSQL') 
    print('%d samples, for brick %s' % (len(Samp),brickname))
    print('First 2 sources have: ')
    for sam in Samp[:2]:
        print('ra=%f, dec=%f' % (sam.ra,sam.dec))
    # Already did these cuts in decals_sim_radeccolors 
    #r0,r1,d0,d1= brickwcs.radec_bounds()
    #Samp.cut( (Samp.ra >= r0)*(Samp.ra <= r1)*\
    #          (Samp.dec >= d0)*(Samp.dec <= d1) )
    # Sort by Sersic n low -> high (if elg or lrg)
    if objtype in ['elg','lrg']:
        if args.cutouts:
            # rhalf ~ 1-2'' at z ~ 1, n~1 
            #Samp=Samp[ (Samp.get('%s_re' % objtype) <= 10.)*\
            #           (Samp.get('%s_n' % objtype) <= 2.) ]
            Samp.set('%s_re' % objtype, np.array([0.5]*len(Samp)))
            Samp.set('%s_n' % objtype, np.array([1.]*len(Samp)))
        else:
            # Usual obiwan
            print('Sorting by sersic n')
            Samp=Samp[np.argsort( Samp.get('%s_n' % objtype) )]
        #    # Dont sort by sersic n for deeplearning cutouts
        #    print('NOT sorting by sersic n')
        #else:
    rowst,rowend= args.rowstart,args.rowstart+maxobjs
    if args.cutouts:
        # Gridded ra,dec for args.stamp_size x stamp_size postage stamps 
        size_arcsec= args.stamp_size * 0.262 * 2 #arcsec, 2 for added buffer
        # 20x20 grid
        dd= size_arcsec / 2. * np.arange(1,21,2).astype(float) #'' offsect from center
        dd= np.concatenate((-dd[::-1],dd))
        dd/= 3600. #arcsec -> deg
        # Don't exceed brick half width - 100''
        assert(dd.max() <= 0.25/2 - 100./3600)
        brickc_ra,brickc_dec= radec_center[0],radec_center[1]
        dec,ra = np.meshgrid(dd+ brickc_dec, dd+ brickc_ra) 
        dec= dec.flatten()
        ra= ra.flatten()
        assert(len(Samp) >= dec.size)
        keep= np.arange(dec.size)
        Samp.cut(keep)
        Samp.set('ra',ra)
        Samp.set('dec',dec)
    # Rowstart -> Rowend
    Samp= Samp[args.rowstart:args.rowstart+maxobjs]
    print('Max sample size=%d, actual sample size=%d' % (maxobjs,len(Samp)))
    assert(len(Samp) <= maxobjs)
    t0= ptime('Got input_sample',t0)

    # Store args in dict for easy func passing
    kwargs=dict(Samp=Samp,\
                brickname=brickname, \
                decals_sim_dir= decals_sim_dir,\
                brickwcs= brickwcs, \
                objtype=objtype,\
                nobj=len(Samp),\
                maxobjs=maxobjs,\
                rowst=rowst,\
                rowend=rowend,\
                do_skipids=args.do_skipids,\
                args=args)

    # Stop if starting row exceeds length of radec,color table
    if len(Samp) == 0:
        fn= get_savedir(**kwargs)+'_exceeded.txt'
        junk= os.system('touch %s' % fn)
        print('Wrote %s' % fn)
        raise ValueError('starting row=%d exceeds number of artificial sources, quit' % rowst)
    
    # Create simulated catalogues and run Tractor
    create_metadata(kwargs=kwargs)
    t0= ptime('create_metadata',t0)
    # do chunks
    #for ith_chunk in chunk_list:
    #log.info('Working on chunk {:02d}/{:02d}'.format(ith_chunk,kwargs['nchunk']-1))
    # Random ra,dec and source properties
    create_ith_simcat(d=kwargs)
    t0= ptime('create_ith_simcat',t0)
    # Run tractor
    do_one_chunk(d=kwargs)
    t0= ptime('do_one_chunk',t0)
    # Clean up output
    if args.cutouts == False:
        do_ith_cleanup(d=kwargs)
        t0= ptime('do_ith_cleanup',t0)
    log.info('All done!')
    return 0
     
if __name__ == '__main__':
    main()
