"""
tests that obiwan runs end to end and get reasonalbe outputs for varietry of
cases. Travis CI runs this script
"""
from __future__ import print_function
if __name__ == "__main__":
    import matplotlib
    matplotlib.use('Agg')
import unittest

import run_200x200_pixel_regions as tools
#from run_200x200_pixel_regions import run_kenobi_main, run_kenobi_main_cosmos
#from run_200x200_pixel_regions import test_flux_shape_measurements
#from run_200x200_pixel_regions import test_detected_simulated_and_real_sources
#from run_200x200_pixel_regions import test_draw_circles_around_sources_check_by_eye

class run_and_analyze(object):
    def __init__(self, survey=None, dataset=None,
                 bands='grz', obj='elg',rowstart=0,
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
        """
        assert(bands in ['z','grz'])
        d= locals()
        del d['self']
        self.params= dict(d)

    def run(self):
        d= dict(self.params)
        if d['checkpoint']:
            # create checkpoint file
            d.update(no_cleanup=True,stage='fitblobs')

        R= tools.run_kenobi_main(**d)
        R.run()

        if d['checkpoint']:
            # restart from checkpoint and finish
            d.update(no_cleanup=False,stage=None)
            R= tools.run_kenobi_main(**d)
            R.run()

    def analyze(self):
        d= dict(self.params)
        if not d['early_coadds']:
            tools.test_flux_shape_measurements(**d)
            tools.test_detected_simulated_and_real_sources(**d)
        tools.test_draw_circles_around_sources_check_by_eye(**d)

class test_main(unittest.TestCase):
    def test_decals(self):
        Test= run_and_analyze(survey='decals',
                              all_blobs=False,on_edge=False,
                              early_coadds=False,
                              checkpoint=False,skip_ccd_cuts=False)

        Test.params.update(dataset='dr5',bands='grz')
        Test.run()
        Test.analyze()

        Test.params.update(bands='z')
        Test.run()
        Test.analyze()

        Test.params.update(early_coadds=True)
        Test.run()
        Test.analyze()

        Test.params.update(early_coadds=False,all_blobs=True)
        Test.run()
        Test.analyze()

        Test.params.update(all_blobs=False,on_edge=True)
        Test.run()
        Test.analyze()

        Test.params.update(on_edge=False,obj='star')
        Test.run()
        Test.analyze()

        # Test.params.update(dataset='dr5',obj='elg',checkpoint=True)
        # Test.run()
        # Test.analyze()

        Test.params.update(dataset='dr3',bands='grz',obj='elg')
        Test.run()
        Test.analyze()

        Test.params.update(bands='z')
        Test.run()
        Test.analyze()

        # Above must simply complete w/out error
        self.assertTrue(True)

    def test_bass_mzls(self):
        Test= run_and_analyze(survey='bass_mzls',obj='elg',
                              skip_ccd_cuts=True,
                              all_blobs=False,on_edge=False,
                              early_coadds=False,
                              checkpoint=False)

        Test.params.update(dataset='dr6',bands='grz')
        Test.run()
        Test.analyze()

        # Above must simply complete w/out error
        self.assertTrue(True)


    def test_cosmos_subset_60(self):
        # runs a 200x200 pixel region but on the full 0.5 GB images
        t= tools.run_kenobi_main_cosmos(survey='decals')
        t.run()
        print('WARNING: no analyze method exists for testcase_cosmos')
        print('adapt run_a_test_case().anlayze()')

        # Above must simply complete w/out error
        self.assertTrue(True)


if __name__ == "__main__":
    unittest.main()
