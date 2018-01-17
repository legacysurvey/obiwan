"""
Continuous integration 'end-to-end' tests of the Obiwan pipeline
"""

from __future__ import print_function
import os
import matplotlib.pyplot as plt
import numpy as np
import re
import sys

import fitsio
import photutils

from obiwan.qa.visual import plotImage, readImage
from obiwan.common import get_brickinfo_hack

try: 
    from astrometry.util.fits import fits_table
    from astrometry.libkd.spherematch import match_radec
    from legacypipe.survey import LegacySurveyData, wcs_for_brick
    from obiwan.kenobi import main,get_parser, get_checkpoint_fn
except ImportError:
    pass

DATASETS= ['dr5','dr3','cosmos']



class Testcase(object):
    """Initialize and run a testcase

    Args:
        name: testcase name
        dataset: string, 'DR3', 'DR5', 
        obj: elg,star
        add_noise: to add Poisson noise to simulated galaxy profiles
        all_blobs: to fit models to all blobs, not just the blobs containing sims
        onedge: to add randoms at edge of region, not well within the boundaries
    """

    def __init__(self, bands='grz',
                 dataset='dr5',
                 obj='elg',rowstart=0,
                 add_noise=False,all_blobs=False,
                 onedge=False, early_coadds=False,
                 checkpoint=False,
                 no_cleanup=False,stage=None):
        assert(dataset in DATASETS)
        self.bands= bands
        self.dataset= dataset
        self.obj= obj
        self.rowstart= rowstart
        self.all_blobs= all_blobs
        self.add_noise= add_noise
        self.onedge= onedge
        self.early_coadds= early_coadds
        self.checkpoint= checkpoint
        self.no_cleanup=no_cleanup
        self.stage=stage

        self.testcase_dir= os.path.join(os.path.dirname(__file__), 
                                        'testcase_%s' % bands)
        self.outname= 'out_testcase_%s_%s_%s' % (bands,
                                self.dataset,self.obj)
        if self.all_blobs:
            self.outname += '_allblobs'
        if self.add_noise:
            self.outname += '_addnoise'
        if self.onedge:
            self.outname += '_onedge'
        if self.early_coadds:
            self.outname += '_coadds'
        self.outdir= os.path.join(os.path.dirname(__file__), 
                                  self.outname)
        self.logfn=os.path.join(self.outdir,'log.txt')
        
        if self.bands == 'grz':
            self.brick='0285m165' 
            self.zoom= [3077, 3277, 2576, 2776]
        elif self.bands == 'z':
            self.brick='1741p242'
            self.zoom= [90, 290, 2773, 2973]
        else:
            raise ValueError('bands= %s no allowed' % bands)

        os.environ["LEGACY_SURVEY_DIR"]= self.testcase_dir
         
    def run(self):
        """run it

        Args:
            no_cleanup: don't run cleanup step
        """
        print('Running testcase: %s' % self.outname)
        extra_cmd_line = []
        if self.add_noise:
            extra_cmd_line += ['--add_sim_noise']
        if self.all_blobs:
            extra_cmd_line += ['--all_blobs']
        if self.early_coadds:
            extra_cmd_line += ['--early_coadds']
        if self.checkpoint:
            extra_cmd_line += ['--checkpoint']
        if self.stage:
            extra_cmd_line += ['--stage',self.stage]
        if self.no_cleanup:
            extra_cmd_line += ['--no_cleanup']

        randoms_fn= os.path.join(os.environ["LEGACY_SURVEY_DIR"], 
                                 'randoms_%s.fits' % self.obj)
        if self.onedge:
            randoms_fn= randoms_fn.replace('.fits','_onedge.fits')


        cmd_line=['--dataset', self.dataset, '-b', self.brick, 
                  '-rs',str(self.rowstart), '-n', '4', 
                  '--zoom', str(self.zoom[0]), str(self.zoom[1]), 
                            str(self.zoom[2]), str(self.zoom[3]),
                  '-o', self.obj, '--outdir', self.outdir,
                  '--randoms_from_fits', randoms_fn] + extra_cmd_line
        parser= get_parser()
        args = parser.parse_args(args=cmd_line)

        if self.checkpoint:
            # Globally redirect stdout
            if os.path.exists(self.logfn):
                os.remove(self.logfn)
            try:
                os.makedirs(os.path.dirname(self.logfn))
            except OSError:
                pass # already exists
            sys.stdout = open(self.logfn, "w")

        main(args=args)

        if self.checkpoint:
            # Reset stdout
            sys.stdout = sys.__stdout__
            # The checkpoint file and log should exist
            ckpt_fn= get_checkpoint_fn(self.outdir,
                              self.brick,self.rowstart)
            assert(os.path.exists(ckpt_fn))
            assert(os.path.exists(self.logfn))




class AnalyzeTestcase(Testcase):
    """Automatically loads the relevant outputs for a given testcase_DR_*

    Args:
        name: like 'testcase_DR5_z_allblobs'

    Attributes:
        brick:
        bands:
        zoom:
        brickwcs:
    """
    def __init__(self, **kwargs):
        super(AnalyzeTestcase, self).__init__(**kwargs)

        # Tolerances
        self.tol= self.get_tolerances()

        survey = LegacySurveyData()
        brickinfo= get_brickinfo_hack(survey,self.brick)
        self.brickwcs = wcs_for_brick(brickinfo)

        self.outdir= os.path.join(os.environ['HOME'],
                           'myrepo/obiwan/tests/end_to_end',
                            self.outname)
        self.rsdir='rs0'

    def get_tolerances(self):
        if self.bands == 'grz':
            mw_trans= 2.e-5 # Not 0 b/c ra,dec of model can vary
            # also amazing agreement
            return {'rhalf':0.65, 
                    'apflux':0.2,
                    'skyflux':2.,
                    'modelflux':4.5,
                    'mw_trans':mw_trans}
        elif self.bands == 'z':
            mw_trans= 5.e-6 # Not 0 b/c ra,dec of model can vary
            if self.add_noise:
                # TODO: tune
                return {'rhalf':0.11, 
                      'apflux':0.25,
                      'skyflux':1.1,
                      'modelflux':6.0,
                      'mw_trans':mw_trans}
            if self.onedge:
                return {'rhalf':0.14, 
                      'apflux':0.2,
                      'skyflux':2.,
                      'modelflux':5.5,
                      'mw_trans':mw_trans}
            if self.obj == 'star':
                # simcat-model amazing agreement
                return {'apflux':0.2,
                        'skyflux':1.1,
                        'modelflux':0.6,
                        'mw_trans':mw_trans}

            
            return {'rhalf':0.11, 
                      'apflux':0.25,
                      'skyflux':1.1,
                      'modelflux':6.0,
                      'mw_trans':mw_trans}


    def load_outputs(self):
        """Each output from the testcase becomes an attribute

        Attributes:
            simcat, obitractor:
            jpg_coadds:
            fits_coadds
        """
        print('Loading from %s' % self.outdir)
        dr= '%s/%s/%s' % (self.brick[:3],self.brick,self.rsdir)
        if not self.early_coadds:
            self.obitractor= fits_table(os.path.join(self.outdir,'tractor',
                                        dr,'tractor-%s.fits' % self.brick))
            self.blobs= fitsio.FITS(os.path.join(self.outdir,'metrics',
                                    dr,'blobs-%s.fits.gz' % self.brick))[0].read()
            self.model_jpg= readImage(os.path.join(self.outdir,'coadd',
                                        dr,'legacysurvey-%s-model.jpg' % self.brick),
                                      jpeg=True)
            self.resid_jpg= readImage(os.path.join(self.outdir,'coadd',
                                            dr,'legacysurvey-%s-resid.jpg' % self.brick),
                                      jpeg=True)
        
        self.simcat= fits_table(os.path.join(self.outdir,'obiwan',
                                    dr,'simcat-%s-%s.fits' % (self.obj,self.brick)))
        
        self.img_jpg= readImage(os.path.join(self.outdir,'coadd',
                                    dr,'legacysurvey-%s-image.jpg' % self.brick),
                                jpeg=True)

        self.img_fits,self.ivar_fits,self.sims_fits= {},{},{}
        for b in self.bands:
            self.img_fits[b]= readImage(os.path.join(self.outdir,'coadd',
                                            dr,'legacysurvey-%s-image-%s.fits.fz' % \
                                             (self.brick,b)))
            self.ivar_fits[b]= readImage(os.path.join(self.outdir,'coadd',
                                            dr,'legacysurvey-%s-invvar-%s.fits.fz' % \
                                              (self.brick,b)))
            self.sims_fits[b]= readImage(os.path.join(self.outdir,'coadd',
                                            dr,'legacysurvey-%s-sims-%s.fits.fz' % \
                                              (self.brick,b)))
            
    def simcat_xy(self):
        """x,y of each simulated source in the fits coadd. Just like the
            bx,by of tractor catalogues
        """
        _,x,y=self.brickwcs.radec2pixelxy(self.simcat.ra,self.simcat.dec)
        self.simcat.set('x',x - self.zoom[0])
        self.simcat.set('y',y - self.zoom[2])

    def match_simcat_tractor(self):
        """matches sim and real sources to tractor cat

        Returns:
            isim,itrac: indices into simcat,tractor 
            ireal: inices into tractor
        """
        # sims to tractor
        rad= 10. * self.brickwcs.pixel_scale() / 3600 #deg
        isim,itrac,d= match_radec(self.simcat.ra, self.simcat.dec, 
                                  self.obitractor.ra, self.obitractor.dec,          
                                  rad,nearest=True)
        # real galaxy to tractor 
        ra_real,dec_real= self.brickwcs.pixelxy2radec(
                                [self.zoom[0] + 100.],
                                [self.zoom[2] + 100.])
        _,ireal,d= match_radec(ra_real, dec_real, 
                               self.obitractor.ra, self.obitractor.dec,          
                               rad,nearest=True)
        return isim,itrac,ireal

    def plots(self):
        """outdir: where write plots to """
        fig,ax=plt.subplots(2,2,figsize=(6,6))
        plotImage().imshow(self.blobs,ax[0,0],qs=None)
        plotImage().imshow(self.img_jpg,ax[0,1],qs=None)
        plotImage().imshow(self.model_jpg,ax[1,0],qs=None)
        plotImage().imshow(self.resid_jpg,ax[1,1],qs=None)
        fn=os.path.join(self.outdir,'blob_img_mod_res.png')
        plt.savefig(fn,dpi=150)
        print('Wrote %s' % fn)
        plt.close()

        fig,ax=plt.subplots(3,3,figsize=(7,6))
        for i,b in enumerate(self.bands) :
            plotImage().imshow(self.img_fits[b],ax[0,i])
            plotImage().imshow(self.sims_fits[b],ax[1,i],qs=None)
            plotImage().circles(self.obitractor.bx,self.obitractor.by,ax[1,i],
                                img_shape=self.sims_fits[b].shape,
                                r_pixels=5./0.262,color='y')
            plotImage().circles(self.simcat.x,self.simcat.y,ax[1,i],
                                img_shape=self.sims_fits[b].shape,
                                r_pixels=4./0.262,color='m')
            plotImage().imshow(self.img_fits[b] - self.sims_fits[b],ax[2,i])
        fn=os.path.join(self.outdir,'fits_coadd_img_sims_res.png')
        plt.savefig(fn,dpi=150)
        print('Wrote %s' % fn)
        plt.close()

    def numeric_tests(self):
        """T: TestcaseOutputs() object """
        isim,itrac,ireal= self.match_simcat_tractor()
        assert((len(isim) == 4) & (len(itrac) == 4))
        if self.all_blobs:
            assert(len(ireal) ==1) # found the real galaxy
        else:
            if self.bands == 'grz':
                assert(len(ireal) == 1) # central galaxy in blob on a sim
            elif self.bands == 'z':
                assert(len(ireal) == 0)

        if not self.obj == 'star':
            rhalf= np.max((self.obitractor[itrac].shapeexp_r,
                           self.obitractor[itrac].shapedev_r),axis=0)
            diff= rhalf - self.simcat.rhalf
            print('delta_rhalf',diff)
            assert(np.all(np.abs(diff) < self.tol['rhalf']))
        # Tractor apflux is nearly bang on to my apflux for sims coadd 
        # plus my apflux for sky in coadd
        # However, Tractor model flux is does not agree with fits coadd counts
        # so its computing on something else and is currently wrong for sim sources
        apers= photutils.CircularAperture((self.simcat.x,self.simcat.y), 
                                           r=3.5/self.brickwcs.pixel_scale())

        for b in self.bands: 
            print('band= %s' % b)

            diff= self.simcat.get('mw_transmission_%s' % b)[isim] -\
                  self.obitractor.get('mw_transmission_%s' % b)[itrac]
            print('delta mw_trans',diff)
            assert(np.all(np.abs(diff) < self.tol['mw_trans'])) 

            apy_table = photutils.aperture_photometry(self.img_fits[b], apers)
            img_apflux= np.array(apy_table['aperture_sum'])
            apy_table = photutils.aperture_photometry(self.sims_fits[b], apers)
            sims_apflux= np.array(apy_table['aperture_sum'])
            obitractor_apflux= self.obitractor[itrac].get('apflux_'+b)[:,5]
            # my apflux vs tractor apflux
            diff= img_apflux - obitractor_apflux
            print('delta_apflux',diff) 
            assert(np.all(np.abs(diff) < self.tol['apflux']))
            # sky flux is small
            diff= img_apflux - sims_apflux
            print('delta_skyflux',diff) 
            assert(np.all(np.abs(diff) < self.tol['skyflux']))
            # tractor model flux within 5-6 nanomags of input flux 
            diff= self.simcat.get(b+'flux') -\
                    self.obitractor[itrac].get('flux_'+b)
            print('delta_modelflux',diff) 
            assert(np.all(np.abs(diff) < self.tol['modelflux']))
        print('passed: Numeric Tests')

    def qualitative_tests(self):
        """T: TestcaseOutputs() object """
        if self.checkpoint:
            # log file is assumed to exist and it must have
            # skipped the correct number of blobs
            assert(self.logfn)
            assert(os.path.exists(self.logfn))
            with open(self.logfn,'r') as foo:
                text= foo.read()
            foundIt= re.search(r'Skipping\s4\sblobs\sfrom\scheckpoint\sfile', text)
            assert(foundIt)
        print('passed: Qualitative Tests')


def test_case(dataset='dr5',
               z=True,grz=False,
               obj='elg',
               add_noise=False,all_blobs=False,
               onedge=False, early_coadds=False,
               checkpoint=False):
    """
    Args:
        dataset: dr5, dr3, cosmos
        z, grz: to run the z and/or grz testcases
        all_blobs: to fit models to all blobs, not just the blobs containing sims
        add_noise: to add Poisson noise to simulated galaxy profiles
        onedge: to add randoms at edge of region, not well within the boundaries
        early_coadds: write coadds before model fitting and stop there
        dataset: no reason to be anything other than DR5 for these tests
    """
    d= dict(obj=obj,dataset=dataset,
            add_noise=add_noise,all_blobs=all_blobs,
            onedge=onedge,early_coadds=early_coadds,
            checkpoint=checkpoint)   
    if z:
        bands= 'z'
    elif grz:
        bands= 'grz'
    d.update(bands=bands) 

    if checkpoint:
        # create checkpoint file
        d.update(no_cleanup=True,stage='fitblobs')

    T= Testcase(**d)
    T.run()
    
    if checkpoint:
        # restart from checkpoint and finish
        d.update(no_cleanup=False,stage=None)
        T= Testcase(**d)
        T.run()
     

    if not early_coadds:
        A= AnalyzeTestcase(**d)
        #if not checkpoint:
        # checkpoint doesn't run cleanup
        A.load_outputs()
        A.simcat_xy()
        A.plots()
        A.numeric_tests()
        A.qualitative_tests()


def test_main():
    """travis CI"""
    d=dict(dataset='dr5',
           z=True,grz=False,
           obj='elg',
           all_blobs=False,onedge=False,
           early_coadds=False,
           checkpoint=False)

    # dr5
    d.update(early_coadds=True)
    test_case(**d)
    
    d.update(early_coadds=False)
    test_case(**d)
    d.update(all_blobs=True)
    test_case(**d)

    d.update(all_blobs=False,onedge=True)
    test_case(**d)

    d.update(obj='star',onedge=False)
    test_case(**d)

    d.update(obj='elg',z=False,grz=True)
    test_case(**d)

    d.update(z=True,grz=False,
             all_blobs=False,
             checkpoint=True)
    test_case(**d)

    # dr3
    d.update(dataset='dr3',
             z=True,grz=False,
             obj='elg',
             all_blobs=False,onedge=False,
             early_coadds=False,
             checkpoint=False)
    test_case(**d)


if __name__ == "__main__":
    #test_dataset_DR3()
    #test_dataset_DR5()
    # Various tests can do 
    d=dict(dataset='dr5',
           z=True,grz=False,
           obj='elg',
           all_blobs=False,onedge=False,
           early_coadds=False,
           checkpoint=False)

    #test_main()
    
    #d.update(early_coadds=True)
    #test_case(**d)

    test_case(**d)
    
