"""
Continuous integration 'end-to-end' tests of the Obiwan pipeline
"""
from __future__ import print_function
import unittest
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
from obiwan.kenobi import main,get_parser, get_checkpoint_fn
# except ImportError:
    # pass

SURVEYS= ['decals','bass_mzls']
DATASETS= ['dr5','dr3','cosmos','dr6']


def nanomag2mag(nmgy):
    return -2.5 * (np.log10(nmgy) - 9)

class run_obiwan_200x200_region(object):
    """Injects 4 sources into a 200x200 pixel image and runs legacypipe

    Works for either DECaLS or MzLS/BASS. Use Analyze_Testcase() on the result

    Args:
        survey: decals or bass_mzls
        name: testcase name
        dataset: string, 'DR3', 'DR5',
        obj: elg,star
        add_noise: to add Poisson noise to simulated galaxy profiles
        all_blobs: to fit models to all blobs, not just the blobs containing sims
        on_edge: to add randoms at edge of region, not well within the boundaries
    """

    def __init__(self, survey=None, dataset=None,
                 bands='grz', obj='elg',rowstart=0,
                 add_noise=False,all_blobs=False,
                 on_edge=False, early_coadds=False,
                 checkpoint=False,
                 no_cleanup=False,stage=None,
                 skip_ccd_cuts=False):
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

        main(args=args)

        if self.checkpoint:
            # Reset stdout
            sys.stdout = sys.__stdout__
            # The checkpoint file and log should exist
            ckpt_fn= get_checkpoint_fn(self.outdir,
                              self.brick,self.rowstart)
            assert(os.path.exists(ckpt_fn))
            assert(os.path.exists(self.logfn))


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
        if obj == 'star':
            mw_trans= 5.e-6
            modelflux=0.6
            return {'modelflux':modelflux,
                    'mw_trans': mw_trans}
        else:
            modelflux=6.0
            if bands == 'z':
                mw_trans= 5.e-6
                rhalf= 0.15
                return {'rhalf':rhalf,
                        'modelflux':modelflux,
                        'mw_trans':mw_trans}
            else:
                mw_trans= 2.e-5
                rhalf= 0.65
                return {'rhalf':rhalf,
                          'modelflux':modelflux,
                          'mw_trans':mw_trans}

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


class setup_tests(run_obiwan_200x200_region):
    """Loads the outputs from run_obiwan_200x200_region().run() and measurement tolerances

    Args:
        kwargs: the same inputs to run_obiwan_200x200_region

    Attributes:
        same as run_obiwan_200x200_region
        tol: tolerance dict for flux, rhalf, etc.
    """
    def __init__(self,**kw):
        super(setup_tests, self).__init__(**kw)
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

class test_flux_shape_measurements(setup_tests):
    def __init__(self,**kw):
        super(test_flux_shape_measurements, self).__init__(**kw)
        self.load_outputs()
        self.run_test()

    def load_outputs(self):
        """Each output from the testcase becomes an attribute

        Attributes:
            simcat, obitractor:
            jpg_coadds:
            fits_coadds
        """
        print('Loading from %s' % self.config_dir)
        self.randoms= fits_table(os.path.join(self.config_dir,'randoms_elg.fits'))
        print('Loading from %s' % self.outdir)
        dr= '%s/%s/%s' % (self.brick[:3],self.brick,self.rsdir)
        self.obitractor= fits_table(os.path.join(self.outdir,'tractor',
                                    dr,'tractor-%s.fits' % self.brick))
        self.simcat= fits_table(os.path.join(self.outdir,'obiwan',
                                    dr,'simcat-%s-%s.fits' % (self.obj,self.brick)))

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


class test_detected_simulated_and_real_sources(setup_tests):
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


class test_draw_circles_around_sources_check_by_eye(setup_tests):
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
        print("test_draw_circles_around_sources_check_by_eye\nLOOK AT THE PNGs")


def test_case(survey=None,dataset=None,
               z=True,grz=False, obj='elg',
               add_noise=False,all_blobs=False,
               on_edge=False, early_coadds=False,
               checkpoint=False, skip_ccd_cuts=False):
    """
    Args:
        survey: one of SURVEYS
        dataset: one of DATASETS
        z, grz: to run the z and/or grz testcases
        all_blobs: to fit models to all blobs, not just the blobs containing sims
        add_noise: to add Poisson noise to simulated galaxy profiles
        on_edge: to add randoms at edge of region, not well within the boundaries
        early_coadds: write coadds before model fitting and stop there
        dataset: no reason to be anything other than DR5 for these tests
    """
    d= dict(survey=survey,dataset=dataset,
            obj=obj,
            add_noise=add_noise,all_blobs=all_blobs,
            on_edge=on_edge,early_coadds=early_coadds,
            checkpoint=checkpoint,
            skip_ccd_cuts=skip_ccd_cuts)
    if z:
        bands= 'z'
    elif grz:
        bands= 'grz'
    d.update(bands=bands)

    if checkpoint:
        # create checkpoint file
        d.update(no_cleanup=True,stage='fitblobs')

    T= run_obiwan_200x200_region(**d)
    T.run()

    if checkpoint:
        # restart from checkpoint and finish
        d.update(no_cleanup=False,stage=None)
        T= run_obiwan_200x200_region(**d)
        T.run()


    if not early_coadds:
        A= Analyze_Testcase(**d)
        #if not checkpoint:
        # checkpoint doesn't run cleanup
        A.load_outputs()
        A.simcat_xy()
        A.plots()
        A.numeric_tests()
        A.qualitative_tests()


def test_main():
    """This is what travis CI runs"""

    # BASS/MzLS
    d=dict(survey='bass_mzls',dataset='dr6',
           z=False,grz=True, skip_ccd_cuts=True,
           obj='elg',
           all_blobs=False,on_edge=False,
           early_coadds=False,
           checkpoint=False)
    test_case(**d)

    # DECaLS
    d.update(survey='decals',dataset='dr5',
             z=True,grz=False,
             skip_ccd_cuts=False)

    # dr5
    d.update(early_coadds=True)
    test_case(**d)

    d.update(early_coadds=False)
    test_case(**d)
    d.update(all_blobs=True)
    test_case(**d)

    d.update(all_blobs=False,on_edge=True)
    test_case(**d)

    d.update(obj='star',on_edge=False)
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
             all_blobs=False,on_edge=False,
             early_coadds=False,
             checkpoint=False)
    test_case(**d)

    d.update(grz=True)
    test_case(**d)


class TestcaseCosmos(object):
    def __init__(self, survey=None,
                 dataset='cosmos',subset=60):
        self.survey=survey
        self.dataset=dataset
        self.subset=subset
        self.rowstart=0
        self.obj='elg'
        if self.survey == 'decals':
            self.brick='1501p020'
        else:
            raise ValueError('survey = bass_mzls not supported yet')
        self.outdir= os.path.join(os.path.dirname(__file__),
                                  'out_testcase_%s_cosmos_subset%d' % \
                                  (survey,self.subset))
        os.environ["LEGACY_SURVEY_DIR"]= os.path.join(os.path.dirname(__file__),
                                                'testcase_cosmos')


    def run(self):
        randoms_fn= os.path.join(os.environ["LEGACY_SURVEY_DIR"],
                                 'randoms_%s.fits' % self.obj)
        cmd_line=['--subset', str(self.subset),
                  '--dataset', self.dataset, '-b', self.brick,
                  '-rs',str(self.rowstart), '-n', '4',
                  '-o', self.obj, '--outdir', self.outdir,
                  '--randoms_from_fits', randoms_fn]
        parser= get_parser()
        args = parser.parse_args(args=cmd_line)

        main(args=args)



if __name__ == "__main__":
    #test_main()
    d=dict(survey='decals',dataset='dr5',
           z=True,grz=False,
           obj='elg',
           all_blobs=False,on_edge=False,
           early_coadds=False,
           checkpoint=False,
           skip_ccd_cuts=False)

    d.update(bands='grz')
    for key in ['grz','z']:
        del d[key]
    #test_flux_truth_vs_measured(**d)

    # test_flux_shape_measurements(**d)
    # test_detected_simulated_and_real_sources(**d)
    # test_draw_circles_around_sources_check_by_eye(**d)
    # raise ValueError('decals done')

    d.update(survey='bass_mzls',dataset='dr6',
             skip_ccd_cuts=True,
             z=False,grz=True)
    #test_case(**d)

    d.update(bands='grz')
    for key in ['grz','z']:
        del d[key]

    test_flux_shape_measurements(**d)
    test_detected_simulated_and_real_sources(**d)
    test_draw_circles_around_sources_check_by_eye(**d)
    raise ValueError('bass_mzls done')
    #test_flux_truth_vs_measured(**d)



    if False:
        d.update(bands='grz')
        for key in ['grz','z']:
            del d[key]
        A= Analyze_Testcase(**d)
        A.load_outputs()
        A.simcat_xy()
        A.plots()
        A.numeric_tests()
        A.qualitative_tests()

    #t= TestcaseCosmos(survey='decals')
    #t.run()
