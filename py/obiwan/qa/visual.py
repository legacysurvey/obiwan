import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle,Wedge
from matplotlib.collections import PatchCollection
import os
import skimage.io
import fitsio

from astrometry.util.fits import fits_table
from legacypipe.survey import LegacySurveyData, wcs_for_brick

class plotImage(object):
    """Helper functions for displaying image and overlaying circles around sources
    
    Args:
        img: need to give as initial input b/c some helper funcs that dont directly
            use img, need its shape at least, see circles()
    """
    def imshow(self,img,ax,qs=[0.5,99.5]):
        if img.shape[-1] == 3:
            #minmax=np.percentile(np.sum(img,axis=2),q=qs)
            minmax=[None,None]
            cmap=None
        elif qs is None:
            minmax=[None,None]
            cmap='gray'
        else:
            minmax=np.percentile(img,q=qs)
            cmap='gray'
        ax.imshow(img, interpolation='none', origin='lower',
                  cmap=cmap,vmin=minmax[0],vmax=minmax[1])
        ax.tick_params(direction='out')
        
    def circles(self,xs,ys,ax,
                img_shape=None,
                xslice=None,yslice=None,
                r_pixels=5./0.262,color='y'):
        """
        xs,ys: x,y positions of sources in pixels, e.g. tractor.bx or simcat.x
        img_shape: needed when xslice or yslice is None
        xlice,yslice: slice() objects into the image array
        r_pixels: radius circles in pixels
        """
        if (xslice is None) | (yslice is None):
            assert(not img_shape is None)
        if xslice is None:
            xslice= slice(0,img_shape[0])
        if yslice is None:
            yslice= slice(0,img_shape[1])
        keep= self.justInSlice(xs,ys,xslice,yslice)
        xpos,ypos= xs[keep]-xslice.start,ys[keep]-yslice.start
        
        dr= r_pixels/ 20 
        patches=[Wedge((x, y), r_pixels + dr, 0, 360,dr)
                 for x,y in zip(xpos, ypos) ]
        coll = PatchCollection(patches, color=color) #,alpha=1)
        ax.add_collection(coll)
        
    def justInSlice(self,x,y,xslice,yslice):
        """Returns bool array of x,y positions in the slice()"""
        return ((x >= xslice.start) & 
               (x <= xslice.stop) &
               (y >= yslice.start) & 
               (y <= yslice.stop))

def readImage(fn,jpeg=False,ext=1):
    """Reads FITS and jpeg images so that x,y indices refer to the same pixels
    regardless of image format. x,y and fits correspond so the jpeg is rotated and flipped 
    to align with fits
    
    Args:
        fn: image filename
        jpeg: bool, is is a jpeg?
    """
    if jpeg:
        img= skimage.io.imread(fn)
        for i in range(3):
            img[:,:,i]= np.rot90(img[:,:,i].T,1)
    else:
        img= fitsio.FITS(fn)[ext].read()
    return img

def sliceImage(img,
               xslice=slice(None,None),yslice=slice(None,None)):
    """Not sure why, but simcat.x[xslice],simcat.y[yslice]
    corresponds to img[yslice,xslice], eg inverted for the image"""
    return img[yslice,xslice,...]

def flux2mag(flux):
    return -2.5*np.log10(1e-9 * flux)



class TestcaseOutputs(object):
    """Automatically loads the relevant outputs for a given testcase_DR_*

    Args:
        name: like 'testcase_DR5_z_allblobs'

    Attributes:
        brick:
        bands:
        zoom:
        brickwcs:
    """
    def __init__(self,name):
        self.name= name
        if '_grz' in name:
            self.brick='0285m165' 
            self.bands= ['g','r','z']
            self.zoom= [3077, 3277, 2576, 2776]
        else:
            self.brick='1741p242'
            self.bands= ['z']
            self.zoom= [90, 290, 2773, 2973]
        os.environ["LEGACY_SURVEY_DIR"]= os.path.join(os.environ['HOME'],
                                            'myrepo/obiwan/tests/end_to_end',
                                            name.replace('_allblobs',''))

        survey = LegacySurveyData()
        brickinfo = survey.get_brick_by_name(self.brick)
        self.brickwcs = wcs_for_brick(brickinfo)
    
    def load(self):
        """Each output from the testcase becomes an attribute

        Attributes:
            simcat, obitractor:
            jpg_coadds:
            fits_coadds
        """
        OUT_DIR= os.path.join(os.environ['HOME'],
                           'myrepo/obiwan/tests/end_to_end',
                           'out_'+self.name,'elg/%s/%s/rs0/' % \
                             (self.brick[:3],self.brick)
                         )
        self.simcat= fits_table(os.path.join(OUT_DIR,'obiwan/simcat-elg-%s.fits' % self.brick))
        self.obitractor= fits_table(os.path.join(OUT_DIR,'tractor/tractor-%s.fits' % self.brick))

        self.blobs= fitsio.FITS(os.path.join(OUT_DIR,'metrics/blobs-%s.fits.gz' % self.brick))[0].read()

        self.img_jpg= readImage(os.path.join(OUT_DIR,'coadd/legacysurvey-%s-image.jpg' % self.brick),
                           jpeg=True)
        self.model_jpg= readImage(os.path.join(OUT_DIR,'coadd/legacysurvey-%s-model.jpg' % self.brick),
                             jpeg=True)
        self.fresid_jpg= readImage(os.path.join(OUT_DIR,'coadd/legacysurvey-%s-resid.jpg' % self.brick),
                             jpeg=True)

        self.img_fits,self.ivar_fits,self.sims_fits= {},{},{}
        for b in self.bands:
            self.img_fits[b]= readImage(os.path.join(OUT_DIR,'coadd/legacysurvey-%s-image-%s.fits.fz' % \
                                             (self.brick,b)))
            self.ivar_fits[b]= readImage(os.path.join(OUT_DIR,'coadd/legacysurvey-%s-invvar-%s.fits.fz' % \
                                              (self.brick,b)))
            self.sims_fits[b]= readImage(os.path.join(OUT_DIR,'coadd/legacysurvey-%s-sims-%s.fits.fz' % \
                                              (self.brick,b)))
            
    def simcat_xy(self):
        """x,y of each simulated source in the fits coadd. Just like the
            bx,by of tractor catalogues
        """
        _,x,y=self.brickwcs.radec2pixelxy(self.simcat.ra,self.simcat.dec)
        self.simcat.set('x',x - self.zoom[0])
        self.simcat.set('y',y - self.zoom[2])

