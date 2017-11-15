import numpy as np
from obiwan.qa.visual import readImage,sliceImag

HDF5_KEYS= ['g','r','z','gr','gz','rz','grz']

class Cutouts(object):
    """Object for exracting all sims and real galaxy cutouts 

        Args:
            outdir: '/home/kaylan/myrepo/obiwan/tests/end_to_end/out_testcase_DR5_grz_allblobs'
            ls_dir: LEGACY_SURVEY_DIR
    """

    def __init__(self,outdir,legacy_survey_dir=None):
        self.outdir= outdir
        if ls_dir is None:
            os.environ["LEGACY_SURVEY_DIR"]= os.path.join(os.path.dirname(__file__), 
                                                          'testcase_DR5_z')
        self.survey = LegacySurveyData()
 
    def get_images_in_bricks(bricks, zoom=None):
        """Write the hdf5 image files for all rs/* in these bricks

        Args:
            bricks: list of bricks
            zoom: if legacypipe was run with zoom option
        """
        for brick in bricks:
            # One hdf5 file for each brick
            self.hdf5_dict={}

            rs_dirs= glob(os.path.join(self.outdir,
                                'elg/%s/%s/*rs0*' % \
                            (brick[:3],brick)
                          )
            assert(len(rs_dirs) > 0))
            self.get_brickwcs(brick)
            for rs_dir in rs_dirs:
                self.load_data(rs_dir,brick)
                self.simcat_xy(zoom=zoom)
                self.extract()
        

    def load_data(self,rs_dir,brick):
        """loads all necessary info for each brick,rs_dir combination

        Args:
            rs_dir: path/to/rs0, rs300, rs300_skipid, etc
        """
        print('Loading from %s' % rs_dir)
        self.simcat= fits_table(os.path.join(rs_dir,'obiwan/simcat-elg-%s.fits' % brick))
        self.obitractor= fits_table(os.path.join(rs_dir,'tractor/tractor-%s.fits' % brick))

        self.img_fits,self.ivar_fits,self.sims_fits= {},{},{}
        self.bands= (pd.Series( glob(os.path.join(rs_dir,'coadd/*-image-*.fits.fz')) )
                     .str.replace('.fits.fz','')
                     .str[-1].values)
        assert(self.bands.size > 0) 
        for b in self.bands: 
            self.img_fits[b]= readImage(os.path.join(rs_dir,'coadd/legacysurvey-%s-image-%s.fits.fz' % \
                                             (brick,b)))
            self.ivar_fits[b]= readImage(os.path.join(rs_dir,'coadd/legacysurvey-%s-invvar-%s.fits.fz' % \
                                              (brick,b)))

    def get_brickwcs(self,brick):
        brickinfo = survey.get_brick_by_name(brick)
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
        key= ''.join(sorted(self.bands)) # like 'rz'
        assert(key in HDF5_KEYS)
        for cat in self.simcat:
            xslc= slice(int(cat.x)-hw,int(cat.x)+hw)
            yslc= slice(int(cat.y)-hw,int(cat.y)+hw)
            # N x N x Number of bands
            self.hdf5_dict[key][str(cat.id)+'/img']= \
                np.array([sliceImg(self.img_fits[band], xslice=xslc,yslice=yslc)
                          for band in key]).T
            self.hdf5_dict[key][str(cat.id)+'/ivar']= \
                np.array([sliceImg(self.ivar_fits[band], xslice=xslc,yslice=yslc)
                          for band in key]).T
                         
    def save_hdf5(self):
        for bands in hdf5.keys():

        self.hdf5= os.path.join(self.outdir,
                                       'elg/%s/%s' % (brick[:3],brick),
                                       ''



