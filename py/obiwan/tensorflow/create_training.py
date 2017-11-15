import numpy as np
import os
from glob import glob
import h5py
import pandas as pd

from astrometry.util.fits import fits_table

from legacypipe.survey import LegacySurveyData, wcs_for_brick
from obiwan.qa.visual import readImage,sliceImage

HDF5_KEYS= ['g','r','z','gr','gz','rz','grz']

class BrickStamps(object):
    """Object for exracting all sims and real galaxy cutouts 

        Args:
            ls_dir: LEGACY_SURVEY_DIR, like 'tests/end_to_end/testcase_DR5_grz'
            obj_dir: '/home/kaylan/myrepo/obiwan/tests/end_to_end/out_testcase_DR5_grz_allblobs'

    """

    def __init__(self,ls_dir=None,obj_dir=None):
        self.obj_dir= obj_dir
        os.environ["LEGACY_SURVEY_DIR"]= ls_dir
        self.survey = LegacySurveyData()
 
    def write_hdf5_for_brick(self,brick, zoom=None):
        """Write the hdf5 image files for all rs/* in this brick

        Args:
            brick:
            zoom: if legacypipe was run with zoom option
        """
        self.get_brickwcs(brick)
        rs_dirs= glob(os.path.join(self.obj_dir,'elg/%s/%s/*rs0*' % \
                                   (brick[:3],brick)))
        assert(len(rs_dirs) > 0)
        # set of bands in this brick
        self.bands= (pd.Series( glob(os.path.join(rs_dirs[0],'coadd/*-image-*.fits.fz')) )
                     .str.replace('.fits.fz','')
                     .str[-1].values)
        assert(self.bands.size > 0)
        # One hdf5 file for this brick
        # like '_rz.hdf5'
        self.bands_str= ''.join(sorted(self.bands))
        assert(self.bands_str in HDF5_KEYS)
        hdf5_fn= os.path.join(self.obj_dir,
                              'elg/%s/%s/' % (brick[:3],brick),
                              'img_ivar_%s.hdf5' % self.bands_str)
        self.hdf5_obj = h5py.File(hdf5_fn, "w")
        # Many rs*/ dirs per brick
        for rs_dir in rs_dirs:
            self.load_brick(rs_dir,brick)
            self.simcat_xy(zoom=zoom)
            self.extract()
        self.hdf5_obj.close()
        
    def load_brick(self,rs_dir,brick):
        """loads all necessary info for each brick,rs_dir combination

        Args:
            rs_dir: path/to/rs0, rs300, rs300_skipid, etc
        """
        print('Loading from %s' % rs_dir)
        self.simcat= fits_table(os.path.join(rs_dir,'obiwan/simcat-elg-%s.fits' % brick))
        self.obitractor= fits_table(os.path.join(rs_dir,'tractor/tractor-%s.fits' % brick))

        self.img_fits,self.ivar_fits,self.sims_fits= {},{},{} 
        for b in self.bands: 
            self.img_fits[b]= readImage(os.path.join(rs_dir,'coadd/legacysurvey-%s-image-%s.fits.fz' % \
                                             (brick,b)))
            self.ivar_fits[b]= readImage(os.path.join(rs_dir,'coadd/legacysurvey-%s-invvar-%s.fits.fz' % \
                                              (brick,b)))

    def get_brickwcs(self,brick):
        brickinfo = self.survey.get_brick_by_name(brick)
        self.brickwcs = wcs_for_brick(brickinfo)

    def simcat_xy(self,zoom=None):
        """x,y of each simulated source in the fits coadd. Just like the
            bx,by of tractor catalogues
        """
        _,x,y=self.brickwcs.radec2pixelxy(self.simcat.ra,self.simcat.dec)
        if zoom:
        	x -= zoom[0]
        	y -= zoom[2]
        self.simcat.set('x',x)
        self.simcat.set('y',y)

    def extract(self,hw=20):
        """

        Args:
            hw: half-width, pixels, (hw*2) x (hw*2) image cutout
        """
        for cat in self.simcat:
            xslc= slice(int(cat.x)-hw,int(cat.x)+hw)
            yslc= slice(int(cat.y)-hw,int(cat.y)+hw)
            # N x N x Number of bands                                                        
            _ = self.hdf5_obj.create_dataset(str(cat.id)+'/img', 
                                             chunks=True, \
                data= np.array([sliceImage(self.img_fits[band],
                                           xslice=xslc,yslice=yslc)
                                for band in self.bands_str]).T)
            _ = self.hdf5_obj.create_dataset(str(cat.id)+'/ivar', 
                                             chunks=True, \
                data= np.array([sliceImage(self.ivar_fits[band],
                                           xslice=xslc,yslice=yslc)
                                for band in self.bands_str]).T)


if __name__ == '__main__':
    name= 'testcase_DR5_grz'
    if '_grz' in name:
        brick='0285m165' 
        zoom= [3077, 3277, 2576, 2776]
    else:
        brick='1741p242'
        zoom= [90, 290, 2773, 2973]

    repo_dir= '/home/kaylan/myrepo/obiwan/'
    ls_dir= os.path.join(repo_dir,
                         'tests/end_to_end',name)
    obj_dir= os.path.join(repo_dir,
                         'tests/end_to_end','out_'+name)
    
    Stamps= BrickStamps(ls_dir=ls_dir, obj_dir=obj_dir)
    for brick in [brick]:
        Stamps.write_hdf5_for_brick(brick, zoom=zoom)

