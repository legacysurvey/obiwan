import numpy as np
import os
from glob import glob
import h5py
import pandas as pd

from astrometry.util.fits import fits_table

from legacypipe.survey import LegacySurveyData, wcs_for_brick
from obiwan.qa.visual import readImage,sliceImage
from obiwan.common import dobash,get_rsdir
from obiwan.runmanager.status import get_final_dir

import galsim

HDF5_KEYS= ['g','r','z','gr','gz','rz','grz']

def flux2mag(nmgy):
    return -2.5 * (np.log10(nmgy) - 9)

class SimStamps(object):
    """Object for exracting sim cutouts 

        Args:
            ls_dir: LEGACY_SURVEY_DIR, like 'tests/end_to_end/testcase_DR5_grz'
            outdir: path to dir containing obiwan,coadd,tractor dirs
    """

    def __init__(self,outdir=None,
                 ls_dir=None):
        """outdir: required
           ls_dir: not needed if env var LEGACY_SURVEY_DIR already set
        """
        self.outdir= outdir
        if ls_dir:
            os.environ["LEGACY_SURVEY_DIR"]= ls_dir
        self.survey = LegacySurveyData()

    def get_brickwcs(self,brick):
        brickinfo = self.survey.get_brick_by_name(brick)
        self.brickwcs = wcs_for_brick(brickinfo)

    def load_data(self,brick,cat_fn,coadd_dir):
        """loads coadd and catalogue data

        Args:
            brick:
            rs_dir: path/to/rs0, rs300, rs300_skipid, etc
        """
        print('Loading from %s' % coadd_dir)
        self.cat= fits_table(cat_fn)

        self.img_fits,self.ivar_fits= {},{}
        for b in self.bands: 
            self.img_fits[b]= readImage(coadd_dir,
                                        'legacysurvey-%s-image-%s.fits.fz' % (brick,b)))
            self.ivar_fits[b]= readImage(coadd_dir,
                                        'legacysurvey-%s-invvar-%s.fits.fz' % (brick,b)))
        # galsim.Image() so can determine overlap w/cutouts
        self.img_gs,self.ivar_gs= {},{}
        for b in self.bands:
            self.img_gs[b]= galsim.Image(self.img_fits[b])
            self.ivar_gs[b]= galsim.Image(self.ivar_fits[b])

    def extract(self,hw=20):
        """For each id,x,y in self.cat, extracts image cutout

        Args:
            hw: half-width, pixels, (hw*2) x (hw*2) image cutout
        """
        for cat in self.cat:
            xslc= slice(cat.x-hw,cat.x+hw)
            yslc= slice(cat.y-hw,cat.y+hw)
            # N x N x Number of bands
            test_img= galsim.Image(np.zeros((2*hw+1,2*hw+1)))
            # y,x because numpy indexing
            test_img.setCenter(cat.x,cat.y)
            olap= test_img.bounds & self.img_gs[self.bands[0]].bounds
            assert(olap.area() > 0)
            if olap.numpyShape() == test_img.array.shape:
                # Grab from fits image b/c aligned better
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

                # _ = self.hdf5_obj.create_dataset(str(cat.id)+'/img', 
                #                                  chunks=True, \
                #     data= np.array([self.img_fits[band][olap].array
                #                     for band in self.bands_str]).T)
                # _ = self.hdf5_obj.create_dataset(str(cat.id)+'/ivar', 
                #                                  chunks=True, \
                #     data= np.array([self.ivar_fits[band][olap].array
                #                     for band in self.bands_str]).T)

            else:
                # On edge
                # Note, galsim.Image() cannot be 3D
                img= [test_img.copy()]*len(self.bands)
                ivar= [test_img.copy()]*len(self.bands)
                for i,band in enumerate(self.bands_str):
                    img[i][olap] += self.img_gs[band][olap]
                    ivar[i][olap] += self.ivar_gs[band][olap]

                _ = self.hdf5_obj_onedge.create_dataset(str(cat.id)+'/img', 
                                                        chunks=True, \
                    data= np.array([d.array for d in img]).T)
                _ = self.hdf5_obj_onedge.create_dataset(str(cat.id)+'/ivar', 
                                                        chunks=True, \
                    data= np.array([d.array for d in ivar]).T)

    def run(self,brick,zoom=None):
        """Write the hdf5 image files for all rs/* in this brick

        Args:
            brick: brickname
            zoom: if legacypipe was run with zoom option
        """
        self.get_brickwcs(brick)
        assert(len(rs_dirs) > 0)
        # coadd fits images must exist
        self.set_paths_to_data(brick)
        coadd_fns= glob(os.path.join(self.coadd_dirs[0],
                                     '*-image-*.fits.fz'))
        assert(len(coadd_fns) > 0)
        # set of bands in this brick
        self.bands= (pd.Series(coadd_fns)
                     .str.replace('.fits.fz','')
                     .str[-1].values)
        assert(self.bands.size > 0)
        # One hdf5 file for this brick, like '_rz.hdf5'
        self.bands_str= ''.join(sorted(self.bands))
        assert(self.bands_str in HDF5_KEYS)
        self.set_hdf5_fns(brick,self.bands_str)
    
        if os.path.exists(self.hdf5_fn):
            print('Skipping %s, hdf5 already exists: %s' % (brick,self.hdf5_fn))
            return None
        self.hdf5_obj = h5py.File(self.hdf5_fn, "w")
        self.hdf5_obj_onedge = h5py.File(self.hdf5_fn_onedge, "w")
        # Many rs*/ dirs per brick
        for cat_fn,coadd_dir in zip(self.cat_fns,self.coadd_dirs):
            self.load_data(brick,cat_fn,coadd_dir)
            self.add_xyid(zoom=zoom)
            self.apply_cuts()
            self.extract()
        self.hdf5_obj.close()
        self.hdf5_obj_onedge.close()
        print('Wrote %s' % hdf5_fn)
        print('Wrote %s'% hdf5_fn_onedge)
        
    def set_paths_to_data(self,brick):
        """lists of catalogues filenames and coadd dirs"""
        rs_dirs= glob(os.path.join(self.outdir,'coadd',
                                   brick[:3],brick,'*rs*'))
        rs_dirs= [os.path.basename(a)
                  for a in rs_dirs]
        for rs_dir in rs_dirs:
            self.cat_fns= os.path.join(self.outdir,'obiwan',
                                       brick[:3],brick,rs_dir,
                                       'simcat.fits')
            self.coadd_dirs= os.path.join(self.outdir,'coadd',
                                       brick[:3],brick,rs_dir)

    def set_hdf5_fns(self,brick,bands_str):
        dr= os.path.join(self.outdir,'hdf5',
                         brick[:3],brick)
        self.hdf5_fn= os.path.join(dr,
                                   'img_ivar_%s.hdf5' % bands_str)
        self.hdf5_fn_onedge= self.hdf5_fn.replace('.hdf5',
                                        '_onedge.hdf5')
        try:
            dobash('mkdir -p %s' % dr)
        except ValueError:
            print('hdf5 dir already exists: ',dr)

    def set_xyid(self,zoom=None):
        _,x,y=self.brickwcs.radec2pixelxy(self.cat.ra,self.cat.dec)
        if zoom:
            x -= zoom[0]
            y -= zoom[2]
        self.cat.set('x',x.astype(int))
        self.cat.set('y',y.astype(int))
        assert('id' in self.cat.get_columns())

    def apply_cuts(self):
        pass

        

class TractorStamps(SimStamps):
    def __init__(self,ls_dir=None,outdir=None,
                 savedir=None):
        """Same as SimStamps but for tractor catalogues

        Args:
            savedir: required for tractor not sims b/c cannot write to dr5 dir
        """
        super(TractorStamps,self).__init__(ls_dir=ls_dir,
                                           outdir=outdir)
        self.savedir= savedir
    
    def set_paths_to_data(self,brick):
        """lists of catalogues filenames and coadd dirs"""
        rs_dirs= 
        for rs_dir in rs_dirs:
            self.cat_fns= [os.path.join(self.outdir,'tractor',
                                        brick[:3],
                                        'tractor-%s.fits' % brick)]
            self.coadd_dirs= [os.path.join(self.outdir,'coadd',
                                           brick[:3],brick)]

    def set_hdf5_fns(self,brick,bands_str):
        dr= os.path.join(self.savedir,'hdf5',
                         brick[:3],brick)
        self.hdf5_fn= os.path.join(dr,
                                   'img_ivar_%s.hdf5' % bands_str)
        self.hdf5_fn_onedge= self.hdf5_fn.replace('.hdf5',
                                        '_onedge.hdf5')
        try:
            dobash('mkdir -p %s' % dr)
        except ValueError:
            print('hdf5 dir already exists: ',dr)

    def set_xyid(self,zoom=None):
        x,y=self.cat.bx,self.cat.by
        if zoom:
            x -= zoom[0]
            y -= zoom[2]
        self.cat.set('x',x.astype(int))
        self.cat.set('y',y.astype(int))
        self.cat.set('id',self.cat.objid)

    def apply_cuts(self):
        len_bef= len(self.cat)
        self.cat.cut((self.cat.brick_primary) &
                     (self.cat.nobs_g >= 1) &
                     (self.cat.nobs_r >= 1) &
                     (self.cat.nobs_z >= 1) & 
                     (self.cat.allmask_g == 0) & 
                     (self.cat.allmask_r == 0) & 
                     (self.cat.allmask_z == 0)) 
        self.cat.cut((flux2mag(self.cat.flux_g) > 20.) &
                     (flux2mag(self.cat.flux_r) > 20.) &
                     (flux2mag(self.cat.flux_z) > 19.))
        print('After cut, have %d/%d' % (len(self.cat),len_bef))


def testcase_main(): 
    name= 'testcase_DR5_grz'
    obj='elg'
    onedge=False
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
                         'tests/end_to_end','out_'+name+'_'+obj)
    if onedge:
        obj_dir += '_onedge'
    
    Sim= SimcatStamps(ls_dir=ls_dir, obj_dir=obj_dir)
    for brick in [brick]:
        Sim.run(brick, zoom=zoom)

def mpi_main(nproc=1,which=None,
             outdir=None,ls_dir=None,savedir=None,
             bricks=[]):
    """

    Args:
        nproc: > 1 for mpi4py
        which: one of ['tractor','sim']
        outdir: path to coadd,tractor dirs
        ls_dir: not needed if legacy_survey_dir env var already set
        savedir: only needed if which = tractor, where to write the hdf5 files
        bricks: list bricks to make hdf5 cutouts from
    """
    assert(which in ['tractor','sim'])
    if nproc > 1:
        from mpi4py.MPI import COMM_WORLD as comm
        bricks= np.array_split(bricks, comm.size)[comm.rank]

    if which == 'sim':
        Obj= SimStamps(outdir=outdir,ls_dir=ls_dir)
    else:
        Obj= TractorStamps(outdir=outdir,ls_dir=ls_dir,
                           savedir=savedir)
    
    for brick in bricks:
        Obj.run(brick)


if __name__ == '__main__':
    #testcase_main()

    

    mpi_main(nproc=1,which='sim',
             outdir='/global/cscratch1/sd/kaylanb/obiwan_out/elg_dr5_coadds',
             savedir='',
             ls_dir=None,
             bricks=bricks)

    mpi_main(nproc=1,which='tractor',
             outdir='/global/project/projectdirs/cosmo/data/legacysurvey/dr5',
             savedir='/global/cscratch1/sd/kaylanb/obiwan_out/dr5_cutouts',
             ls_dir=None,
             bricks=bricks)

