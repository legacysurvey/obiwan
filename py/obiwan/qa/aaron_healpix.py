import numpy as np
import os
import healpy as hp
import fitsio
import matplotlib.pyplot as plt
from scipy.stats import sigmaclip
from collections import defaultdict

from astropy.coordinates import Galactic,ICRS
from astropy import units

#from theValidator.catalogues import Matcher,CatalogueFuncs
from astrometry.libkd.spherematch import match_radec
from astrometry.util.fits import fits_table, merge_tables

from obiwan.fetch import fetch_targz
from obiwan.kenobi import dobash

DOWNLOAD_ROOT = "http://portal.nersc.gov/project/desi/users/kburleigh/"

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
    def __init__(self,targz_dir, allDR5=False):
        self.targz_dir= targz_dir
        if allDR5:
          self.drname='dr5_qa_70k'
        else:
          self.drname='dr5_qa'

    def fetch(self):
        name='healpix.tar.gz'
        fetch_targz(os.path.join(DOWNLOAD_ROOT,self.drname,name),
                    os.path.joni(self.targz_dir,self.drname))
        
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
        hdu=fitsio.FITS(os.path.join(self.targz_dir,self.drname,fn))
        data=hdu[0].read()
        nmatch= hdu[1].read()
        return data,nmatch

    def get_radec(self,data,keep):
        """Return ra,dec,subset_data healpix arrays for boolean array keep"""
        nside= Healpix().get_nside( len(data) )
        ra,dec= hp.pix2ang(nside,np.where(keep)[0],lonlat=True)
        return ra,dec, data[keep]

class EmptyClass(object):
    pass

class footprint_wSFD(object):
  """makes nice figure showing DECaLS,MzLS,BASS footprint ontop of sfd98 dust
  
  Example: 
    Foot= footprint_wSFD('/home/kaylan/mydata')
    sfd= Foot.get_footprint_object()
    Foot.plot_footprint_object()
  """
  def __init__(self,data_dir='/home/kaylan/mydata'):
    self.data_dir= data_dir
    self.map_dir= os.path.join(data_dir,'sfd98')
    self.tile_dir= os.path.join(data_dir,'svn_tiles')
    # Download data
    self.download_sfd98_healpix()
    self.download_decals_mzls_tiles()

  def get_footprint_object(self):
    """Returns footprint object 'sfd'"""
    # work with SFD map and Decals/Mzls tiles
    # lonlat from SFD healpix is in galactic coords, convert this to Celestial
    hdu=fitsio.FITS(os.path.join(self.map_dir,'lambda_sfd_ebv.fits'))
    sfd= EmptyClass()
    temp= hdu[1].read()
    sfd.temp= temp['TEMPERATURE'] 
    npix= Healpix().get_nside(len(sfd.temp))
    assert(npix == 512)
    sfd.l_indeg,sfd.b_indeg= hp.pix2ang(512,np.where(sfd.temp > 0)[0],nest=True,lonlat=True)
    #inPlane= np.where((sfd_gal_dec > -20) & (sfd_gal_dec < 20))[0]
    trans= Galactic(l=sfd.l_indeg * units.degree,
                    b=sfd.b_indeg * units.degree)
    radec= trans.transform_to(ICRS)
    sfd.ra,sfd.dec= radec.ra.value, radec.dec.value
    
    all_tiles=fits_table(os.path.join(self.tile_dir,'mosaic-tiles_obstatus.fits'))
    wdes_tiles=fits_table(os.path.join(self.tile_dir,'decam-tiles_obstatus.fits'))
    
    inDESI= ( (all_tiles.in_desi_orig == 1) |
              (all_tiles.in_desi == 1))
    inDecals= ( (inDESI) &
                (all_tiles.dec <= 30.)) 
                #(mzls_decals.in_des == 0))
    inMzls=   ( (inDESI) &
                (all_tiles.dec > 30.)) 
                #(mzls_decals.in_des == 0))
    inDes= (  (wdes_tiles.in_desi_orig == 1) |
              (wdes_tiles.in_desi == 1))
    inDes=    ( (inDes) &
                (wdes_tiles.in_des == 1))
    #above30= mzls.dec > 30.
    #inDESI= ( (mzls.in_desi_orig == 1) |
    #          (mzls.in_desi == 1))
    #inMzls= ( (above30) &
    #          (inDESI))

    #desi= merge_tables([mzls,decals],columns='fillzero')
    des= wdes_tiles.copy()
    del wdes_tiles
    des.cut(inDes)
    mzls= all_tiles.copy()
    decals= all_tiles.copy()
    del all_tiles
    mzls.cut(inMzls)
    decals.cut(inDecals)
    
    ps= Healpix().get_pixscale(len(sfd.temp),unit='deg')
    # match_radec(ref,obs): for each point in ref, return matching point in obs
    print('matching tiles to healpix centers')
    I,J,d= match_radec(mzls.ra,mzls.dec, sfd.ra,sfd.dec, ps*8)
    sfd.ipix_mzls= list( set(J) )

    I,J,d= match_radec(decals.ra,decals.dec, sfd.ra,sfd.dec, ps*8)
    sfd.ipix_decals= list( set(J) )

    I,J,d= match_radec(des.ra,des.dec, sfd.ra,sfd.dec, ps*8)
    sfd.ipix_des= list( set(J) )
    
    return sfd

    # legasurvey pts fill in in ps*3
    I,J,d= match_radec(legsurvey.ra,legsurvey.dec, sfd.ra,sfd.dec, ps*3)
    sfd.ipix_legsurvey= set(J)
    # des fills in with ps*8
    I,J,d= match_radec(des.ra,des.dec, sfd.ra,sfd.dec, ps*8)
    sfd.ipix_legsurvey.union( set(J) )
    sfd.ipix_legsurvey= list(sfd.ipix_legsurvey)
    return sfd

  def plot_footprint_object(self,footprint_obj):
      """sfd is what footprint_wSFD() returns"""
      temp= np.log10(footprint_obj.temp)
      temp[footprint_obj.ipix_legsurvey]= 2.
      hp.mollview(temp,nest=True,flip='geo',title='Mollweide Projection, Galactic Coordinates',unit='',max=-0.5)
      hp.graticule(c='k',lw=1) 
      plt.savefig('footprint_wSFD.png',dpi=150)

  def modify_healpy_colorbar1():
      pass
  
  def modify_healpy_colorbar2():
      x, y, z = np.random.random((3, 30))
      z = z * 20 + 0.1
      
      # Set some values in z to 0...
      z[:5] = 0
      
      cmap = plt.get_cmap('jet', 20)
      cmap.set_under('gray')
      
      fig, ax = plt.subplots()
      cax = ax.scatter(x, y, c=z, s=100, cmap=cmap, vmin=0.1, vmax=z.max())
      fig.colorbar(cax, extend='min')

  def download_sfd98_healpix(self):
    """downloads data if isnt on computer"""
    tar_name= 'sfd98.tar.gz'
    map_name= 'sfd98/lambda_sfd_ebv.fits'
    if not os.path.exists(self.map_dir):
      os.makedirs(self.map_dir)
      fetch_targz(os.path.join(DOWNLOAD_ROOT,'obiwan',tar_name),
                  self.data_dir)

  def download_decals_mzls_tiles(self):
    """downloads data if isnt on computer"""
    tar_name= 'svn_tiles.tar.gz'
    mosaic_nm= 'mosaic-tiles_obstatus.fits'
    decals_nm= 'decam-tiles_obstatus.fits'
    if not os.path.exists(self.tile_dir):
      os.makedirs(self.tile_dir)
      fetch_targz(os.path.join(DOWNLOAD_ROOT,'obiwan',tar_name),
                  self.data_dir)


class Bricks(object):
    def __init__(self,targz_dir,decals=True):
        self.bricks= fits_table(os.path.join(targz_dir,
                                'legacysurveydir','survey-bricks.fits.gz'))
        if decals:
            self.bricks.cut( (self.bricks.dec > -30) &
                             (self.bricks.dec < 30))

    def get_nearest_brick(self,ra,dec):
        """given an ra,dec returns the nearest brick"""
        deg_per_brick=0.25
        #imatch,imiss,d2d= Matcher().match_within(heal,brick, dist= deg_per_brick/2)
        I,J,d = match_radec(ra,dec, self.bricks.ra,self.bricks.dec,
                            deg_per_brick, nearest=True)
        return self.bricks.brickname[ J[0] ]

    def get_nearest_bricks(self,ra_list,dec_list):
        bricks=[]
        for ra,dec in zip(ra_list,dec_list):
            bricks.append( self.get_nearest_brick(ra,dec) )
        return bricks


class Plots(object):
    def __init__(self,outdir='./', close=True):
        self.outdir= outdir
        self.close= close
    
    def basic(self,data,min=None,max=None):
        hp.mollview(data,min=min,max=max,nest=False)
        
    def mollzoom(self, ra,dec,hp_vals,name,
                 vlim=None,ralim=None,declim=None,
                 figsize=(5,5)):
        plt.figure(figsize=figsize)
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
            #                            (bri,brick,brickv))
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
  parser.add_argument('--allDR5', action="store_true", default=False,help='only on NERSC, load GB sized healpix maps with data from all of DR5')
  args = parser.parse_args()
  print(args)
  
  d= Data(args.targz_dir, allDR5=args.allDR5)
  d.fetch()
  data,nmatch= d.get_data(psf_or_aper=args.psf_or_aper,
                          which=args.which)
  #nside= Healpix().get_nside( len(data) )
  #deg_per_healpix= get_pixscale(npix,unit='deg')
  
  p= Plots(close=False,outdir=args.outdir)
  ra,dec,cut_data= d.get_radec(data, keep= data != 0)
  p.mollzoom(ra,dec,np.abs(cut_data),'test',
             ralim=[ra.min(),ra.max()],declim=[dec.min(),dec.max()],
             vlim=(0.01,0.1))

  # Find closest brick to each of 5 largest deviations
  hi_to_low= np.sort(np.abs(cut_data))[::-1]
  worst= defaultdict(list)
  for data_val in hi_to_low[:10]:
      i= np.where( np.abs(cut_data) == data_val)[0]
      worst['dev'].append( data_val )
      worst['ra'].append( ra[i][0] )
      worst['dec'].append( dec[i][0] )

  B= Bricks(args.targz_dir)
  worst['brick']= B.get_nearest_bricks(worst['ra'],worst['dec'])
  with open('worst_%s_%s.txt' % (args.psf_or_aper,args.which),'w') as f:
      for dev,ra,dec,brick in zip(worst['dev'],
                                  worst['ra'],worst['dec'],worst['brick']):
          f.write('%.2f %.2f %.2f %s\n' % (dev,ra,dec,brick))
        
  
  
  #orig_code(data,nmatch)
  #b=fits_table("/global/cscratch1/sd/kaylanb/dr5_qa/brick_table_psfhealpixdecam-ps1-0128-ddec.fits")
  #get_DR5_ccds(b.brickname)    
  #raise ValueError



