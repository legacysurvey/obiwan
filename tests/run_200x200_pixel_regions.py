"""
Similar to legacypipe/py/test/runbrick_test.py, this runs and analyzes the
output of legacypipe running on 200x200 pixels regions that obiwan
has injected simulated sources intoself.

See 'test_200x200_pixel_regions.py' for how to use this file
"""
from __future__ import print_function
if __name__ == "__main__":
    import matplotlib
    matplotlib.use('Agg')
import os
import matplotlib.pyplot as plt
import numpy as np
import re
import sys

import fitsio

from obiwan.qa.visual import plotImage, readImage
from obiwan.common import get_brickinfo_hack

# try:
from astrometry.util.fits import fits_table
from astrometry.libkd.spherematch import match_radec
from legacypipe.survey import LegacySurveyData, wcs_for_brick
from obiwan.kenobi import get_parser, get_checkpoint_fn
from obiwan.kenobi import main as main_kenobi
# except ImportError:
    # pass

SURVEYS= ['decals','bass_mzls']
DATASETS= ['dr3','dr5','dr6','cosmos']


def nanomag2mag(nmgy):
    return -2.5 * (np.log10(nmgy) - 9)

class run_kenobi_main(object):
    """Runs main() in kenobi.py, which runs legacypipe

    Works for either DECaLS or MzLS/BASS. Use setup_testcase() on the result

    Args:
        survey: decals or bass_mzls
        dataset: string, 'DR3', 'DR5',
        bands: 'grz','z'
        obj: 'elg','star'
        rowstart: default is 0 because reads randoms table from 1st row
        add_noise: to add Poisson noise to simulated galaxy profiles
        all_blobs: to fit models to all blobs, not just the blobs containing sims
        on_edge: to add randoms at edge of region, not well within the boundaries
        early_coadds: creates coadds and stops. useful for ML training/test samples
        checkpoint: whether to save model fitting (fitblobs) checkpoints
        skip_ccd_cuts= True to use all ccds in the survey-ccds fits tables
        no_cleanup: at the end of running, obiwan reorganzes outputs, True means don't do this
        stage: runbrick.py stages, default is None which runs everything
    """

    def __init__(self, survey=None, dataset=None,
                 bands='grz', obj='elg',rowstart=0,
                 add_noise=False,all_blobs=False,
                 on_edge=False, early_coadds=False,
                 checkpoint=False,skip_ccd_cuts=False,
                 no_cleanup=False,stage=None):
        assert(survey in SURVEYS)
        assert(dataset in DATASETS)
        self.survey= survey
        self.dataset= dataset
        self.bands= bands
        self.obj= obj
        self.rowstart= rowstart
        self.all_blobs= all_blobs
        self.add_noise= add_noise
        self.on_edge= on_edge
        self.early_coadds= early_coadds
        self.checkpoint= checkpoint
        self.no_cleanup=no_cleanup
        self.stage=stage
        self.skip_ccd_cuts= skip_ccd_cuts

        self.testcase_dir= os.path.join(os.path.dirname(__file__),
                                        'testcase_%s_%s' % (survey,bands))
        self.outname= 'out_testcase_%s_%s_%s_%s' % (survey,bands,
                                self.dataset,self.obj)
        if self.all_blobs:
            self.outname += '_allblobs'
        if self.add_noise:
            self.outname += '_addnoise'
        if self.on_edge:
            self.outname += '_on_edge'
        if self.early_coadds:
            self.outname += '_coadds'
        if self.checkpoint:
            self.outname += '_ckpt'
        self.outdir= os.path.join(os.path.dirname(__file__),
                                  self.outname)
        self.logfn=os.path.join(self.outdir,'log.txt')

        if self.survey == 'decals' and self.bands == 'grz':
            self.brick='0285m165'
            self.zoom= [3077, 3277, 2576, 2776]
        elif self.survey == 'decals' and self.bands == 'z':
            self.brick='1741p242'
            self.zoom= [90, 290, 2773, 2973]
        elif self.survey == 'bass_mzls':
            self.brick='2176p330'
            self.zoom= [2776,2976, 2900,3100]
        else:
            raise ValueError('bands= %s no allowed' % bands)

        os.environ["LEGACY_SURVEY_DIR"]= self.testcase_dir

    def run(self):
        """Add the sources and run legacypipe

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
        if self.skip_ccd_cuts:
            extra_cmd_line += ['--skip_ccd_cuts']

        randoms_fn= os.path.join(os.environ["LEGACY_SURVEY_DIR"],
                                 'randoms_%s.fits' % self.obj)
        if self.on_edge:
            randoms_fn= randoms_fn.replace('.fits','_on_edge.fits')


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

        main_kenobi(args=args)

        if self.checkpoint:
            # Reset stdout
            sys.stdout = sys.__stdout__
            # The checkpoint file and log should exist
            ckpt_fn= get_checkpoint_fn(self.outdir,
                              self.brick,self.rowstart)
            assert(os.path.exists(ckpt_fn))
            assert(os.path.exists(self.logfn))

class run_kenobi_main_cosmos(object):
    """Supports subsets 60-69"""
    def __init__(self, survey=None):
        self.dataset='cosmos'
        self.survey=survey
        self.rowstart=0
        self.obj='elg'
        if self.survey == 'decals':
            self.brick='1501p020'
        else:
            raise ValueError('survey = bass_mzls not supported yet')
        os.environ["LEGACY_SURVEY_DIR"]= os.path.join(os.path.dirname(__file__),
                                                'testcase_cosmos')


    def run(self):
        self.subset=60
        print('WARNING: testcase_cosmos runs subset 60, only, for simplicity')
        print('but obiwan works for subsets 60-69 when running on a full dataset')
        self.outdir= os.path.join(os.path.dirname(__file__),
                                  'out_testcase_%s_cosmos_subset%d' % \
                                  (self.survey,self.subset))

        # download data
        if not os.path.exists(os.environ["LEGACY_SURVEY_DIR"]):
            print('Do the following to download the testcase, they are 500 MB')
            print('cd %s' % os.path.dirname(__file__))
            print('wget http://portal.nersc.gov/project/desi/users/kburleigh/obiwan/testcase_cosmos.tar.gz')
            print('tar -xzvf testcase_cosmos.tar.gz')
            print("testcase_cosmos is also on kaylanb's hpss at nersc")
            return 0 #sys.exit(0) # exit for success
        randoms_fn= os.path.join(os.environ["LEGACY_SURVEY_DIR"],
                                 'randoms_%s.fits' % self.obj)
        cmd_line=['--subset', str(self.subset),
                  '--dataset', self.dataset, '-b', self.brick,
                  '-rs',str(self.rowstart), '-n', '4',
                  '-o', self.obj, '--outdir', self.outdir,
                  '--randoms_from_fits', randoms_fn]
        parser= get_parser()
        args = parser.parse_args(args=cmd_line)

        main_kenobi(args=args)

class Tolerances(object):
    @staticmethod
    def get(survey=None,obj=None,bands=None):
        assert survey in SURVEYS
        assert(bands in ['grz','z'])
        assert(obj in ['elg','star'])
        return getattr(Tolerances(),survey)(obj,bands)
        #elif survey == 'bass_mzls':
        #    return Tolerances().bass_mzls(obj,bands)


    @staticmethod
    def decals(obj=None,bands=None):
        """Returns dict of Tolerances

        Args:
            bands: either 'grz' or 'z'
            obj: either 'elg' or 'star'
            add_noise: was Poisson noise added to the simulated source profile?
            on_edge: did the test case inject simulated sources onto CCD edges?
        """
        tol= dict(mw_transmission_input_minus_measured=2.e-5,
                  mag_psql_minus_input_corrected_for_ext= 1e-14,
                  mag_input_minus_measured=dict(g=0.12,
                                                r=0.1,
                                                z=0.2),
                  rhalf_input_minus_measured= 0.15)
        if bands == 'grz':
            tol.update(rhalf_input_minus_measured=0.35)
        if obj == 'star':
            tol.update(mag_input_minus_measured=dict(g=0.02,
                                                     r=0.02,
                                                     z=0.02))
        return tol
        # if bands == 'z':
        #     rhalf= 0.15
        # else:
        #     rhalf= 0.65

    @staticmethod
    def bass_mzls(obj=None,bands=None):
        """Returns dict of Tolerances

        Args:
            bands: 'grz'
            obj: either 'elg' or 'star'
        """
        return dict(mw_transmission_input_minus_measured=5.e-6,
                    mag_psql_minus_input_corrected_for_ext= 1e-14,
                    mag_input_minus_measured=dict(g=2.6,
                                                  r=2.6,
                                                  z=0.11),
                    rhalf_input_minus_measured= 0.1)


class analysis_setup(run_kenobi_main):
    """Loads the outputs from run_kenobi_main().run() and measurement tolerances

    Args:
        kwargs: the same inputs to run_kenobi_main

    Attributes:
        same as run_kenobi_main
        tol: tolerance dict for flux, rhalf, etc.
    """
    def __init__(self,**kw):
        super(analysis_setup, self).__init__(**kw)
        # raise ValueError('hey')
        self.tol= Tolerances().get(survey=self.survey,
                                   bands=self.bands,obj=self.obj)
        self.config_dir= os.path.join(os.path.dirname(self.outdir),
                                      'testcase_%s_%s' % \
                                        (kw['survey'],kw['bands']))
        self.rsdir='rs0'

        survey = LegacySurveyData()
        brickinfo= get_brickinfo_hack(survey,self.brick)
        self.brickwcs = wcs_for_brick(brickinfo)

    def load_outputs(self):
        """Loads every output we could possibly need to evaluate a test

        Attributes:
            simcat, obitractor:
            jpg_coadds:
            fits_coadds
        """
        print('Loading from %s' % self.config_dir)
        self.randoms= fits_table(os.path.join(self.config_dir,'randoms_elg.fits'))
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

def at_most_N_is_false(bool_array,N=1):
    """Returns True if bool_array has at most N elements that are False

    Example, all(bool_array) is equivalent to at_most_N_is_false(bool_array,N=0)

    Args:
        N: the number of false elements allowed"""
    return bool_array.sum() >= len(bool_array)-1

class test_flux_shape_measurements(analysis_setup):
    def __init__(self,**kw):
        super(test_flux_shape_measurements, self).__init__(**kw)
        self.load_outputs()
        self.run_test()

    # def load_outputs(self):
    #     """Each output from the testcase becomes an attribute
    #
    #     Attributes:
    #         simcat, obitractor:
    #         jpg_coadds:
    #         fits_coadds
    #     """
    #     print('Loading from %s' % self.outdir)
    #     dr= '%s/%s/%s' % (self.brick[:3],self.brick,self.rsdir)
    #     self.obitractor= fits_table(os.path.join(self.outdir,'tractor',
    #                                 dr,'tractor-%s.fits' % self.brick))
    #     self.simcat= fits_table(os.path.join(self.outdir,'obiwan',
    #                                 dr,'simcat-%s-%s.fits' % (self.obj,self.brick)))

    def run_test(self):
        """Compare input flux and shape parameters to Tractor's"""
        isim,itrac,d = match_radec(self.simcat.ra,self.simcat.dec,
                                   self.obitractor.ra,self.obitractor.dec,
                                   1./3600,nearest=True)

        print('mw_transmission_input_minus_measured')
        for b in self.bands:
            diff= self.simcat.get('mw_transmission_%s' % b)[isim] -\
                  self.obitractor.get('mw_transmission_%s' % b)[itrac]
            print(b+' :',diff)
            assert(np.all(np.abs(diff) < self.tol['mw_transmission_input_minus_measured']))
        print('PASS')

        print('mag_psql_minus_input_corrected_for_ext')
        for b in self.bands:
            dmag= self.randoms.get(b) - \
                    nanomag2mag(self.simcat.get(b+'flux')/self.simcat.get('mw_transmission_'+b))
            print(b+' :',dmag)
            if (self.on_edge) or (self.obj == 'star'):
                assert(np.all(np.abs(dmag) < 1.e-3))
            else:
                assert(np.all(np.abs(dmag) < self.tol['mag_psql_minus_input_corrected_for_ext']))
        print('PASS')

        print('mag_input_minus_measured')
        for b in self.bands:
            dmag= nanomag2mag(self.simcat.get(b+'flux')[isim]) - \
                    nanomag2mag(self.obitractor.get('flux_'+b)[itrac])
            print(b+': ',dmag)
            assert(np.all(np.abs(dmag) < self.tol['mag_input_minus_measured'][b]))
        print('PASS')

        print('rhalf_input_minus_measured')
        if self.obj != 'star':
            rhalf= np.max((self.obitractor[itrac].shapeexp_r,
                           self.obitractor[itrac].shapedev_r),axis=0)
            diff= self.simcat.rhalf - rhalf
            print(diff)
            assert(at_most_N_is_false(np.abs(diff) < self.tol['rhalf_input_minus_measured'],
                                      N=1))
        print('PASS')


class test_detected_simulated_and_real_sources(analysis_setup):
    def __init__(self,**kw):
        super(test_detected_simulated_and_real_sources, self).__init__(**kw)
        self.load_outputs()
        self.run_test()

    def get_index_of_real_galaxy_at_center(self):
        """There is a real galaxy at center, which tractor cat index is it?"""
        # real galaxy to tractor
        ra_real,dec_real= self.brickwcs.pixelxy2radec(
                                [self.zoom[0] + 100.],
                                [self.zoom[2] + 100.])
        _,ireal,d= match_radec(ra_real, dec_real,
                               self.obitractor.ra, self.obitractor.dec,
                               10./3600,nearest=True)
        return ireal

    def run_test(self):
        isim,itrac,d = match_radec(self.simcat.ra,self.simcat.dec,
                            self.obitractor.ra,self.obitractor.dec,
                            1./3600,nearest=True)
        print('all_input_sources_recovered_by_tractor')
        assert((len(isim) == 4) & (len(itrac) == 4))
        print('PASS')

        ireal= self.get_index_of_real_galaxy_at_center()
        print('real_galaxy_at_center_recovered_by_tractor')
        if self.all_blobs:
            assert(len(ireal) ==1) # found the real galaxy
        else:
            if self.bands == 'grz':
                assert(len(ireal) == 1) # central galaxy in blob on a sim
            elif self.bands == 'z':
                assert(len(ireal) == 0)
        print('PASS')


class test_draw_circles_around_sources_check_by_eye(analysis_setup):
    def __init__(self,**kw):
        super(test_draw_circles_around_sources_check_by_eye, self).__init__(**kw)
        self.load_outputs()
        # Add x,y pix to simcat table, juts like tractor has bx,by
        _,x,y=self.brickwcs.radec2pixelxy(self.simcat.ra,self.simcat.dec)
        self.simcat.set('x',x - self.zoom[0])
        self.simcat.set('y',y - self.zoom[2])
        # make the plots
        self.run_test()

    def run_test(self):
        if not self.early_coadds:
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
            print("test_draw_circles_around_sources_check_by_eye\nLOOK AT THE PNGs")
        else:
            for name in ['tractor','metrics']:
                assert(not os.path.exists(os.path.join(self.outdir,name)))
            for name in ['obiwan','coadd']:
                assert(os.path.exists(os.path.join(self.outdir,name)))
            dr= '%s/%s/%s' % (self.brick[:3],self.brick,self.rsdir)
            for name in ['model','resid']:
                assert(not os.path.exists(os.path.join(self.outdir,'coadd',
                            dr,'legacysurvey-%s-%s.jpg' % (self.brick,name))))
            for name in ['image']:
                assert(os.path.exists(os.path.join(self.outdir,'coadd',
                            dr,'legacysurvey-%s-%s.jpg' % (self.brick,name))))
            print("test_draw_circles_around_sources_check_by_eye: PASS")
