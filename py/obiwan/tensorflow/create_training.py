import numpy as np

class Cutouts(object):
    """Initialize and run a testcase

    Args:
        name: testcase name
        dataset: string, 'DR3', 'DR5', etc
        zoom: the zoom array pass to runbrick
        all_blobs,add_noise: booleans
    """
    def __init__(self, brickname,
    			 name='testcase_DR5_z',zoom=None):
        os.environ["LEGACY_SURVEY_DIR"]= os.path.join(os.path.dirname(__file__), 
                                                      name)
        self.brick= brickname
        self.zoom= zoom
        
        survey = LegacySurveyData()
        brickinfo = survey.get_brick_by_name(self.brick)
        self.brickwcs = wcs_for_brick(brickinfo)

    def read(self,outdir='/home/kaylan/myrepo/obiwan/tests/end_to_end/out_testcase_DR5_grz_allblobs'):
        """Each output from the testcase becomes an attribute

        Attributes:
            simcat, obitractor:
            jpg_coadds:
            fits_coadds
        """
        self.outdir= os.path.join(outdir,
        				'elg/%s/%s/rs0/' % \
                        (self.brick[:3],self.brick)
                        )
        print('Loading from %s' % self.outdir)
        self.simcat= fits_table(os.path.join(self.outdir,'obiwan/simcat-elg-%s.fits' % self.brick))
        self.obitractor= fits_table(os.path.join(self.outdir,'tractor/tractor-%s.fits' % self.brick))

        self.img_fits,self.ivar_fits,self.sims_fits= {},{},{}
        for b in self.bands:
            self.img_fits[b]= readImage(os.path.join(self.outdir,'coadd/legacysurvey-%s-image-%s.fits.fz' % \
                                             (self.brick,b)))
            self.ivar_fits[b]= readImage(os.path.join(self.outdir,'coadd/legacysurvey-%s-invvar-%s.fits.fz' % \
                                              (self.brick,b)))
            
    def simcat_xy(self):
        """x,y of each simulated source in the fits coadd. Just like the
            bx,by of tractor catalogues
        """
        _,x,y=self.brickwcs.radec2pixelxy(self.simcat.ra,self.simcat.dec)
        if self.zoom:
        	x -= self.zoom[0]
        	y -= self.zoom[2]
        self.simcat.set('x',x)
        self.simcat.set('y',y)

