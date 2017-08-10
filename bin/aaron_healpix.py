import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from scipy.stats import sigmaclip


def get_nside(num_pix):
    assert(num_pix % 12 == 0)
    return np.sqrt(num_pix / 12)

def get_pixscale(num_pix,unit='deg'):
    assert(unit in ['deg','arcsec'])
    deg2= 4*np.pi * (180/np.pi)**2 / num_pix
    if unit == 'deg':
        return np.sqrt(deg2)
    else:
        return np.sqrt(deg2*3600)

def plot(data,min=None,max=None):
    hp.mollview(data,min=min,max=max,nest=False)

hdu=fitsio.FITS('decam-ps1-0128-ddec.fits')
data=hdu[0].read()
npix= len(data)
nside= get_nside(npix)

_, lo, hi = sigmaclip(data[data != 0], low=3, high=3)
flag= np.logical_or(data < lo, data > hi)
ra,dec= hp.pix2ang(nside,np.where(flag)[0],lonlat=True)

brick= fits_table('survey-bricks.fits.gz')
i_ref,i_obs,dist= match(brick.ra,brick.dec,ra,dec, dist= 0.25/3600)
brick.cut(i_ref)

for brick in bricks:
    T=fits_table()
    ccds= fits_table('%s/coadd/ccds.fits' % brick)
    ccds.set('brick',np.array([brick]*len(ccds)))
    T=merge_tables([T,ccd],fill=zero)

