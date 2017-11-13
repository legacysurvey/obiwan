# tests datasets DR3, DR5, DR3_eBOSS

from __future__ import print_function
import os
import fitsio
import matplotlib.pyplot as plt
import photutils

from astrometry.libkd.spherematch import match_radec

from obiwan.kenobi import main,get_parser
from obiwan.qa.visual import plotImage, readImage
from obiwan.qa.visual import TestcaseOutputs

DATASETS= ['DR3','DR5','DR3_eBOSS']

def plots_for_testcase(outdir,T):
    """T: TestcaseOutputs() object """
    fig,ax=plt.subplots(2,2,figsize=(6,6))
    plotImage().imshow(T.blobs,ax[0,0],qs=None)
    plotImage().imshow(T.img_jpg,ax[0,1],qs=None)
    plotImage().imshow(T.model_jpg,ax[1,0],qs=None)
    plotImage().imshow(T.resid_jpg,ax[1,1],qs=None)
    fn=os.path.join(outdir,'blob_img_mod_res.png')
    plt.savefig(fn,dpi=150)
    print('Wrote %s' % fn)
    plt.close()

    fig,ax=plt.subplots(3,3,figsize=(7,6))
    for i,b in enumerate(T.bands) :
        plotImage().imshow(T.img_fits[b],ax[0,i])
        plotImage().imshow(T.sims_fits[b],ax[1,i],qs=None)
        plotImage().circles(T.obitractor.bx,T.obitractor.by,ax[1,i],
                            img_shape=T.sims_fits[b].shape,
                            r_pixels=5./0.262,color='y')
        plotImage().circles(T.simcat.x,T.simcat.y,ax[1,i],
                            img_shape=T.sims_fits[b].shape,
                            r_pixels=4./0.262,color='m')
        plotImage().imshow(T.img_fits[b] - T.sims_fits[b],ax[2,i])
    fn=os.path.join(outdir,'fits_coadd_img_sims_res.png')
    plt.savefig(fn,dpi=150)
    print('Wrote %s' % fn)
    plt.close()

def numeric_tests(T):
    """T: TestcaseOutputs() object """
    # match sims to tractor and central galaxy to tractor
    rad= 10. * t.brickwcs.pixel_scale() / 3600 #deg
    isim,itrac,d= match_radec(t.simcat.ra, t.simcat.dec, t.obitractor.ra, t.obitractor.dec,          
                              rad,nearest=True)

    # real source at center
    x_real,y_real= 100.,100.
    dpix= np.sqrt((t.obitractor.bx - x_real)**2 + (t.obitractor.by - y_real)**2)       
    ireal= [np.argmin(dpix)]
    assert((len(isim) == 4) & (len(itrac) == 4) & (len(ireal) ==1))

    # tractor cat apflux giving correct flux in coadd fits, but model flux
    # is computing based on something else and currently wrong
    assert()
    # Also, tractor apflux is nearly bang on to my apflux for sims coadd 
    # plus my apflux for sky in coadd, which is what eventually have when
    # tractor model fluxes are correct
    apers= photutils.CircularAperture((T.simcat.x,T.simcat.y), 
                                       r=3.5/T.brickwcs.pixel_scale())
    for b in T.bands: 
        apy_table = photutils.aperture_photometry(T.img_fits[b], apers)
        img_apflux= np.array(apy_table['aperture_sum'])
        apy_table = photutils.aperture_photometry(T.sims_fits[b], apers)
        sims_apflux= np.array(apy_table['aperture_sum'])
        sky_apflux= img_apflux - 
        # measured flux should be sims + sky
        np.all(np.abs(sims_apflux + 
                      sky_apflux - T.obitractor[itrac].get('apflux_'+b)[:,5]) 
                    < 0.5) #nanomags
        # ditto for flux in simcat table
        np.all(np.abs(T.simcat.zflux  + 
                      sky_apflux - T.obitractor[itrac].get('apflux_'+b)[:,5]) 
                    < 1.) #nanomags


def run_dataset(dataset):
    """run a single dataset's test problem

    Args:
        dataset: string, 'DR3', 'DR5', etc
    """
    assert(dataset in DATASETS)
    print('testing dataset: %s' % dataset)
    obj='elg'
    brick='1238p245'
    os.environ["LEGACY_SURVEY_DIR"]= os.path.join(os.path.dirname(__file__), 
                                                  'legacypipedir_%s_dataset_%s' % (brick,dataset))
    outdir = os.path.join(os.path.dirname(__file__),
                          'outdir_%s_%s' % (brick,dataset))
    randoms_from_fits= os.path.join(os.path.dirname(__file__), 
                                    'randoms/%s' % obj,'randoms_%s_10.fits' % brick)

    cmd_line=['--dataset', dataset, '-b', brick, '-n', '2', 
              '-o', 'elg', '--outdir', outdir, '--add_sim_noise',
              '--randoms_from_fits', randoms_from_fits]
    parser= get_parser()
    args = parser.parse_args(args=cmd_line)

    main(args=args)
    assert(True)

def run_testcase(name='testcase_DR5_z',dataset='DR5',
                 zoom=[90, 290, 2773, 2973],
                 all_blobs=False,add_noise=True):
    """run a single dataset's test problem

    Args:
        dataset: string, 'DR3', 'DR5', etc
    """
    assert(dataset in DATASETS)
    print('testcase: %s' % name)
    obj='elg'
    if '_grz' in name:
        brick='0285m165' 
    else:
        brick='1741p242'
    os.environ["LEGACY_SURVEY_DIR"]= os.path.join(os.path.dirname(__file__), 
                                                  name)
    outname= name
    if all_blobs:
        outname += '_allblobs' 

    extra_cmd_line = []
    if add_noise:
        extra_cmd_line += ['--add_sim_noise']
    if all_blobs:
        extra_cmd_line += ['--all_blobs']

    outdir = os.path.join(os.path.dirname(__file__),
                          outname)
    randoms_from_fits= os.path.join(os.path.dirname(__file__), 
                                    name,'randoms.fits')

    cmd_line=['--dataset', dataset, '-b', brick, '-n', '4', 
              '--zoom', str(zoom[0]), str(zoom[1]), str(zoom[2]), str(zoom[3]),
              '-o', 'elg', '--outdir', outdir,
              '--randoms_from_fits', randoms_from_fits] + extra_cmd_line
    parser= get_parser()
    args = parser.parse_args(args=cmd_line)

    main(args=args)

    # Load outputs
    t= TestcaseOutputs(outname)
    t.load()
    t.simcat_xy()

    # Make plots
    # idir='elg/%s/%s/rs0' % (brick[:3],brick)
    # blobs= fitsio.FITS(os.path.join(outdir,idir,
    #                     'metrics/blobs-%s.fits.gz' % brick))[0].read()
    # img_jpg= readImage(os.path.join(outdir,idir,
    #                     'coadd/legacysurvey-%s-image.jpg' % brick),
    #                    jpeg=True)
    # model_jpg= readImage(os.path.join(outdir,idir,
    #                     'coadd/legacysurvey-%s-model.jpg' % brick),
    #                      jpeg=True)
    # resid_jpg= readImage(os.path.join(outdir,idir,
    #                     'coadd/legacysurvey-%s-resid.jpg' % brick),
    #                      jpeg=True)
    plots_for_testcase(outdir,t)

    # Quantitative
    numeric_tests(t)


#fn = os.path.join(outdir, 'tractor', '110', 'tractor-1102p240.fits')
#assert(os.path.exists(fn))
#T = fits_table(fn)
#assert(len(T) == 2)
#print('Types:', T.type)
#assert(T.type[0] == 'REX ')
#cmd = ('(cd %s && sha256sum -c %s)' %
#       (outdir, os.path.join('tractor', '110', 'brick-1102p240.sha256sum')))
#print(cmd)
#rtn = os.system(cmd)
#assert(rtn == 0)

def test_dataset_DR3():
    run_dataset('DR3')
    assert(True)

def test_dataset_DR5():
    run_dataset('DR5')
    assert(True)

def test_cases():
    # Two: z, grz, Possibly modifications: can add sim galaxies at diff 
    #   seaprations (currnetly 4 sep by 12''), Or can can run w/ or wout/ 
    #   all_blobs, etc
    d= dict(name='testcase_DR5_z',dataset='DR5',
            zoom=[90, 290, 2773, 2973])
    #run_testcase(add_noise=False,**d)
    run_testcase(add_noise=False,all_blobs=True,**d)
    
    #d= dict(name='testcase_DR5_grz',dataset='DR5',
    #        zoom=[3077, 3277, 2576, 2776])
    #run_testcase(add_noise=False,**d)
    #run_testcase(add_noise=False,all_blobs=True,**d)
    assert(True)



if __name__ == "__main__":
    #test_dataset_DR3()
    #test_dataset_DR5() 
    test_cases()

    
