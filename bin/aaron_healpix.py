import numpy as np
import os
#import healpy as hp
import fitsio
#import matplotlib.pyplot as plt
from scipy.stats import sigmaclip

from theValidator.catalogues import Matcher,CatalogueFuncs
from astrometry.util.fits import fits_table, merge_tables
from obiwan.fetch import fetch_targz

DOWNLOAD_ROOT = "http://portal.nersc.gov/project/desi/users/kburleigh/dr5_qa/"
OUTDIR= '/home/kaylan/mydata/dr5_qa'

def get_nside(num_pix):
    assert(num_pix % 12 == 0)
    return int( np.sqrt(num_pix / 12) )

def get_pixscale(num_pix,unit='deg'):
    assert(unit in ['deg','arcsec'])
    deg2= 4*np.pi * (180/np.pi)**2 / num_pix
    if unit == 'deg':
        return np.sqrt(deg2)
    else:
        return np.sqrt(deg2*3600)

def plot(data,min=None,max=None):
    hp.mollview(data,min=min,max=max,nest=False)


def get_DR5_ccds(bricknames):
    path='/global/cscratch1/sd/desiproc/DR5_out/'
    T=[]
    for brick in bricknames:
        bri=brick[:3]
        ccd_fn= os.path.join(path,
                             'coadd/%s/%s/legacysurvey-%s-ccds.fits' %
                             (bri,brick,brick))
        try: 
            t=fits_table(ccd_fn)
            t.set('brickname', np.array([brick]*len(t)))
            T.append(t)
            #ccd_fns.append(os.path.join(path,
            #                            'coadd/%s/%s/legacysurvey-%s-ccds.fits' %
            #                            (bri,brick,brick))
        except IOError:
            print('not found: %s' % ccd_fn)
    TT= merge_tables(T, columns='fillzero')
    del T
    savefn= 'brick_allccds.fits'
    TT.writeto(savefn)
    print('Wrote %s' % savefn)
    #CatalogueFuncs().stack(ccd_fns, textfile=False)

# Plots
# hp.mollview(data,min=lo,max=hi)
# data2= data.copy()
# data2[flag == False]=0
# hp.mollview(data2,min=lo,max=hi)
# plt.scatter(brick.ra,brick.dec)

#for brick in bricks:
#    T=fits_table()
#    ccds= fits_table('%s/coadd/ccds.fits' % brick)
#    ccds.set('brick',np.array([brick]*len(ccds)))
#    T=merge_tables([T,ccd],fill=zero)


if __name__ == '__main__':
    b=fits_table("/global/cscratch1/sd/kaylanb/dr5_qa/brick_table_psfhealpixdecam-ps1-0128-ddec.fits")
    get_DR5_ccds(b.brickname)    
    raise ValueError
    name='healpix.tar.gz'
    fetch_targz(DOWNLOAD_ROOT+name,
                os.path.join(OUTDIR,name))
    
    fn='psf/healpix/decam-ps1-0128-ddec.fits'
    hdu=fitsio.FITS(os.path.join(OUTDIR,fn))
    data=hdu[0].read()
    npix= len(data)
    nside= get_nside(npix)
    
    _, lo, hi = sigmaclip(data[data != 0], low=3, high=3)
    flag= np.logical_or(data < lo, data > hi)
    #raise ValueError
    ra,dec= hp.pix2ang(nside,np.where(flag)[0],lonlat=True)
    heal= fits_table()
    for col,arr in zip(['ra','dec'],[ra,dec]):
        heal.set(col, arr)
    
    brick= fits_table('/home/kaylan/mydata/survey-bricks.fits.gz')
    brick.cut( (brick.dec > -40)*\
               (brick.dec < 40))
    #deg_per_healpix= get_pixscale(npix,unit='deg')
    deg_per_brick=0.25
    imatch,imiss,d2d= Matcher().match_within(heal,brick, dist= deg_per_brick/2)
    brick.cut(imatch['obs'])
    id= fn.replace('/','').replace('.fits','')
    savenm= os.path.join(OUTDIR,'brick_table_%s.fits' % id)
    brick.writeto(savenm)
    print('Wrote %s' % savenm)

