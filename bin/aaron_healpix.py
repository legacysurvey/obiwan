import numpy as np
import os
import healpy as hp
import fitsio
import matplotlib.pyplot as plt
from scipy.stats import sigmaclip

#from theValidator.catalogues import Matcher,CatalogueFuncs
from astrometry.libkd.spherematch import match_radec
from astrometry.util.fits import fits_table, merge_tables
from obiwan.fetch import fetch_targz

DOWNLOAD_ROOT = "http://portal.nersc.gov/project/desi/users/kburleigh/dr5_qa/"

class Healpix(object):
    def get_nside(self,num_pix):
        assert(num_pix % 12 == 0)
        return int( np.sqrt(num_pix / 12) )
    
    def get_pixscale(self,num_pix,unit='deg'):
        assert(unit in ['deg','arcsec'])
        deg2= 4*np.pi * (180/np.pi)**2 / num_pix
        if unit == 'deg':
            return np.sqrt(deg2)
        else:
            return np.sqrt(deg2*3600)

class Data(object):
    def __init__(self,targz_dir):
        self.targz_dir= targz_dir
        
    def fetch(self):
        name='healpix.tar.gz'
        fetch_targz(DOWNLOAD_ROOT+name,
                    os.path.join(self.targz_dir,'dr5_qa',name))
        
    def get_data(self,psf_or_aper,which):
        """read healpix data, RING ordered'

        Args:
          psf_or_aper: choices ['psf','aper']
          which: choices ['ddec','dra','g','r','z']

        Returns:
          data: healpix array for data
          nmatch: healpix array for number ps1 matches in each pixel
        """
        fn='%s/healpix/decam-ps1-0128-%s.fits' % (psf_or_aper,which)
        hdu=fitsio.FITS(os.path.join(self.targz_dir,'dr5_qa',fn))
        data=hdu[0].read()
        nmatch= hdu[1].read()
        return data,nmatch

    def get_radec(self,data,keep):
        """Return ra,dec,subset_data healpix arrays for boolean array keep"""
        nside= Healpix().get_nside( len(data) )
        ra,dec= hp.pix2ang(nside,np.where(keep)[0],lonlat=True)
        return ra,dec, data[keep]

class Bricks(object):
    def __init__(self,targz_dir,decals=True):
        self.bricks= fits_table(os.path.join(targz_dir,
                                'legacysurveydir','survey-bricks.fits.gz'))
        if decals:
            self.brick.cut( (brick.dec > -20) & (brick.dec < 30))

    def get_nearest_brick(self,ra,dec):
        """returns nearest brick to given ra,dec"""
        deg_per_brick=0.25
        #imatch,imiss,d2d= Matcher().match_within(heal,brick, dist= deg_per_brick/2)
        I,J,d = match_radec(ra,dec, self.bricks.ra,self.bricks.dec,
                            deg_per_brick/2, nearest=True)
        return self.bricks.brickname[I]


class Plots(object):
    def __init__(self,outdir='./', close=True):
        self.outdir= outdir
        self.close= close
    
    def basic(self,data,min=None,max=None):
        hp.mollview(data,min=min,max=max,nest=False)
        
    def mollzoom(self, ra,dec,hp_vals,name,
                 vlim=None,ralim=None,declim=None):
        plt.figure(figsize=(10,4))
        plt.scatter(ra,dec,c=hp_vals,cmap='rainbow',alpha=0.75)
        if vlim:
            plt.clim(vlim)
        plt.xlabel('Ra');plt.ylabel('Dec')
        if ralim:
            plt.xlim(ralim)
        if declim:
            plt.ylim(declim)
        plt.axes().set_aspect('equal')
        plt.colorbar(orientation='vertical')
        plt.tight_layout()
        fn=os.path.join(self.outdir,name+'.png')
        plt.savefig(fn,dpi=150)
        print('Wrote %s' % fn)
        if self.close:
            plt.close()
            
    def scatter(self, ra,dec,name,
                ralim=None,declim=None):
        plt.figure(figsize=(10,4))
        plt.scatter(ra,dec,c='b',alpha=0.75)
        plt.xlabel('Ra');plt.ylabel('Dec')
        if ralim:
            plt.xlim(ralim)
        if declim:
            plt.ylim(declim)
        plt.axes().set_aspect('equal')
        plt.tight_layout()
        fn=os.path.join(self.outdir,name+'.png')
        plt.savefig(fn,dpi=150)
        print('Wrote %s' % fn)
        if self.close:
            plt.close()

def orig_code(data,nmatch):
    nside= Healpix().get_nside( len(data) )
    
    _, lo, hi = sigmaclip(data[data != 0], low=3, high=3)
    flag= np.logical_or(data < lo, data > hi)
    flag*= (nmatch > 20)
    ra,dec= hp.pix2ang(nside,np.where(flag)[0],lonlat=True)
    
    # PLOTTING
    ralim=[ra.min(),ra.max()]
    declim=[dec.min(),dec.max()]
    my_mollzoom(ra,dec,data[flag],'outliers',
                ralim=ralim,declim=declim, vlim=(lo,hi))
    temp_ra,temp_dec= hp.pix2ang(nside,np.where(np.ones(len(data),bool))[0],lonlat=True)
    keep= (temp_ra >= ralim[0])*\
          (temp_ra <= ralim[1])*\
          (temp_dec >= declim[0])*\
          (temp_dec <= declim[1])
    my_mollzoom(temp_ra[keep],temp_dec[keep],data[keep],'all',
                ralim=ralim,declim=declim, vlim=(lo,hi))
    keep*= (nmatch > 20)
    my_mollzoom(temp_ra[keep],temp_dec[keep],data[keep],'nmatch_gt20',
                ralim=ralim,declim=declim, vlim=(lo,hi))
    
    # Match bricks
    #heal= fits_table()
    #for col,arr in zip(['ra','dec'],[ra,dec]):
    #    heal.set(col, arr)
    brick= fits_table(os.path.join(args.targz_dir,
                                   'legacysurveydir','survey-bricks.fits.gz'))
    brick.cut( (brick.dec > -40)*\
               (brick.dec < 40))
    #deg_per_healpix= get_pixscale(npix,unit='deg')
    deg_per_brick=0.25
    #imatch,imiss,d2d= Matcher().match_within(heal,brick, dist= deg_per_brick/2)
    I,J,d = match_radec(ra,dec, brick.ra,brick.dec,deg_per_brick/2,
                        nearest=True)
    #raise ValueError
    brick.cut(imatch['obs'])
    my_scatter(brick.ra,brick.dec,'bricks',
               ralim=ralim,declim=declim)
    
    id= fn.replace('/','').replace('.fits','')
    savenm= os.path.join(args.outdir,'brick_table_%s.fits' % id)
    brick.writeto(savenm)
    print('Wrote %s' % savenm)    

    
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

#for brick in bricks:
#    T=fits_table()
#    ccds= fits_table('%s/coadd/ccds.fits' % brick)
#    ccds.set('brick',np.array([brick]*len(ccds)))
#    T=merge_tables([T,ccd],fill=zero) 

if __name__ == '__main__':
  from argparse import ArgumentParser
  parser = ArgumentParser()
  parser.add_argument('--targz_dir', type=str, default='/home/kaylan/mydata/dr5_qa',help='where wget data to')
  parser.add_argument('--outdir', type=str, default=None, help='where save plots and fits tables to')
  parser.add_argument('--psf_or_aper', type=str, choices=['psf','aper'],default='psf')
  parser.add_argument('--which', type=str, choices=['ddec','dra','g','r','z'],default='ddec')
  args = parser.parse_args()
  print(args)
  
  d= Data(args.targz_dir)
  d.fetch()
  data,nmatch= d.get_data(psf_or_aper=args.psf_or_aper,
                          which=args.which)
  #nside= Healpix().get_nside( len(data) )
  #deg_per_healpix= get_pixscale(npix,unit='deg')
  
  p= Plots(close=False,outdir=args.outdir)
  ra,dec,cut_data= d.get_radec(data, keep= data != 0)
  p.mollzoom(ra,dec,np.abs(cut_data),'test',ralim=[0,20],declim=[-20,20],
             vlim=(0.01,0.1))

  # Find closest brick to each of 5 largest deviations
  hi_to_low= np.sort(np.abs(cut_data))[::-1]
  for data_val in hi_to_low[:5]:
      i= np.where( np.abs(cut_data) == data_val)[0]
      print(ra[i],dec[i], B.get_nearest_brick(ra[i],dec[i]))
        
  B= Bricks(args.targz_dir)
  print( B.get_nearest_brick(hot_ra,hot_dec) )
  
  
  #orig_code(data,nmatch)
  #b=fits_table("/global/cscratch1/sd/kaylanb/dr5_qa/brick_table_psfhealpixdecam-ps1-0128-ddec.fits")
  #get_DR5_ccds(b.brickname)    
  #raise ValueError



