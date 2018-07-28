"""
saves 64x64 pixel cutouts of each source in a Data Release as HDF5 files
"""

import numpy as np
import os
from glob import glob
import h5py
import pandas as pd

from obiwan.qa.visual import readImage,sliceImage
from obiwan.common import dobash

try:
    from astrometry.util.fits import fits_table
    from legacypipe.survey import LegacySurveyData, wcs_for_brick
    import galsim
except ImportError:
    pass

HDF5_KEYS= ['g','r','z','gr','gz','rz','grz']

def flux2mag(nmgy):
    return -2.5 * (np.log10(nmgy) - 9)

class SimStamps(object):
    """Object for exracting sim cutouts

        Args:
            ls_dir: LEGACY_SURVEY_DIR, like 'tests/end_to_end/testcase_DR5_grz'
            outdir: path to dir containing obiwan,coadd,tractor dirs
    """

    def __init__(self,ls_dir=None,outdir=None,
                 savedir=None, jpeg=False):
        """outdir: required
           ls_dir: not needed if env var LEGACY_SURVEY_DIR already set
           save_dir: where write hdf5 files, outdir if None
        """
        self.outdir= outdir
        self.jpeg= jpeg
        if ls_dir:
            os.environ["LEGACY_SURVEY_DIR"]= ls_dir
        self.savedir= savedir
        if self.savedir is None:
            self.savedir= self.outdir
        self.survey = LegacySurveyData()

    def get_brickwcs(self,brick):
        brickinfo = self.survey.get_brick_by_name(brick)
        self.brickwcs = wcs_for_brick(brickinfo)

    def load_data(self,brick,cat_fn,coadd_dir):
        """loads coadd and catalogue data

        Args:
            brick:
            coadd_dir: path/to/rs0, rs300, rs300_skipid, etc
        """
        print('Loading ra,dec from %s' % (os.path.dirname))
        self.cat= fits_table(cat_fn)

        print('Loading from %s' % coadd_dir)
        self.img_fits,self.ivar_fits= {},{}
        for b in self.bands:
            self.img_fits[b]= readImage(os.path.join(coadd_dir,
                                        'legacysurvey-%s-image-%s.fits.fz' % (brick,b)))
            self.ivar_fits[b]= readImage(os.path.join(coadd_dir,
                                        'legacysurvey-%s-invvar-%s.fits.fz' % (brick,b)))
        self.img_jpeg= readImage(os.path.join(coadd_dir,
                                 'legacysurvey-%s-image.jpg' % (brick)),
                                 jpeg=True)
        # galsim.Image() so can determine overlap w/cutouts
        self.img_gs,self.ivar_gs,self.jpeg_gs= {},{},{}
        for ib,b in enumerate(self.bands):
            self.img_gs[b]= galsim.Image(self.img_fits[b])
            self.ivar_gs[b]= galsim.Image(self.ivar_fits[b])

    def extract(self,hw=32):
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

                _ = self.hdf5_jpeg.create_dataset(str(cat.id)+'/img',
                                                  chunks=True,dtype=np.uint8, \
                    data= sliceImage(self.img_jpeg,
                                     xslice=xslc,yslice=yslc))
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

    def run(self,brick,stampSize=64,applyCuts=True,zoom=None):
        """Write the hdf5 image files for all rs/* in this brick

        Args:
            brick: brickname
            stampSize: height and width in pixes of training image
            zoom: if legacypipe was run with zoom option
        """
        self.get_brickwcs(brick)
        # coadd fits images must exist
        self.set_paths_to_data(brick)
        coadd_fns= glob(os.path.join(self.coadd_dirs[0],
                                     '*-image-*.fits.fz'))
        if len(coadd_fns) == 0:
            raise IOError('no image.fits.fz file here: %s' % self.coadd_dirs[0])
        # set of bands in this brick
        self.bands= (pd.Series(coadd_fns)
                     .str.replace('.fits.fz','')
                     .str[-1].values)
        assert(self.bands.size > 0)
        # One hdf5 file for this brick, like '_rz.hdf5'
        self.bands_str= ''.join(sorted(self.bands))
        assert(self.bands_str in HDF5_KEYS)
        self.set_output_fns(brick,self.bands_str)

        if os.path.exists(self.hdf5_fn):
            try:
                f= h5py.File(self.hdf5_fn, 'r')
                f2= h5py.File(self.hdf5_fn_onedge, 'r')
                f3= h5py.File(self.hdf5_fn_jpeg, 'r')
                if (len(f.keys()) == 0) & (len(f2.keys()) == 0):
                    # remove empty files then make them
                    for fn in [self.hdf5_fn,self.hdf5_fn_onedge,self.hdf5_fn_jpeg]:
                        dobash('rm %s' % fn)
                else:
                    # processing done, skip this brick
                    print('Skipping %s, hdf5 already filled: %s' % (brick,self.hdf5_fn))
                    return None
            except OSError:
                # One of these got messed up, redo it
                for fn in [self.hdf5_fn,self.hdf5_fn_onedge,self.hdf5_fn_jpeg]:
                    os.remove(fn)
                    print('removed ',fn)
        self.hdf5_obj = h5py.File(self.hdf5_fn, "w")
        self.hdf5_obj_onedge = h5py.File(self.hdf5_fn_onedge, "w")
        self.hdf5_jpeg = h5py.File(self.hdf5_fn_jpeg, "w")
        # Many rs*/ dirs per brick
        for cat_fn,coadd_dir in zip(self.cat_fns,self.coadd_dirs):
            if (self.jpeg) & (os.path.basename(coadd_dir) != 'rs0'):
                print('jpeg=True, so skipping %s' % coadd_dir)
                continue
            self.load_data(brick,cat_fn,coadd_dir)
            self.set_xyid(zoom=zoom)
            if applyCuts:
                self.apply_cuts()
                self.write_mag_sorted_ids()
            self.extract(hw=int(stampSize/2))
        self.hdf5_obj.close()
        self.hdf5_obj_onedge.close()
        self.hdf5_jpeg.close()
        print('Wrote %s' % self.hdf5_fn)
        print('Wrote %s' % self.hdf5_fn_onedge)
        print('Wrote %s' % self.hdf5_fn_jpeg)

    def set_paths_to_data(self,brick):
        """lists of catalogues filenames and coadd dirs"""
        search= os.path.join(self.outdir,'coadd',
                             brick[:3],brick,'*rs*',
                             'legacysurvey-%s-ccds.fits' % brick)
        rs_dirs= glob(search)
        rs_dirs= [os.path.basename(os.path.dirname(a))
                  for a in rs_dirs]
        if len(rs_dirs) == 0:
            raise IOError('No rs dirs here: %s' % search)
        self.cat_fns,self.coadd_dirs= \
            zip(*[(os.path.join(self.outdir,'obiwan',
                                brick[:3],brick,rs_dir,
                                'simcat-elg-%s.fits' % brick),
                   os.path.join(self.outdir,'coadd',
                                brick[:3],brick,rs_dir)
                   )
                   for rs_dir in rs_dirs])

    def set_output_fns(self,brick,bands_str):
        dr= os.path.join(self.savedir,'hdf5',
                         brick[:3],brick)
        # hdf5
        self.hdf5_fn= os.path.join(dr,
                                   'img_ivar_%s.hdf5' % bands_str)
        self.hdf5_fn_onedge= self.hdf5_fn.replace('.hdf5',
                                        '_onedge.hdf5')
        self.hdf5_fn_jpeg= self.hdf5_fn.replace('img_ivar_','jpeg_')
        # table of mag sorted ids
        self.sorted_ids_fn= self.hdf5_fn.replace(os.path.basename(self.hdf5_fn),
                                                'sorted_ids.fits')
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
        len_bef=len(self.cat)
        print('After cut, have %d/%d' % (len(self.cat),len_bef))

    def write_mag_sorted_ids(self,band='g'):
        mag= self.get_mag(band)
        inds= np.argsort(mag) # small to large (mag, so brightest to faintest)
        T= fits_table()
        T.set('id',self.cat.id)
        T.set('mag_'+band,mag)
        T= T[inds]
        T.writeto(self.sorted_ids_fn)
        print('Wrote %s' % self.sorted_ids_fn)

    def get_mag(self,band='g'):
        return flux2mag(self.cat.get(band+'flux'))


#######
# Funcs to apply simulated source cuts to tractor catalogues sources

def get_xy_pad(slope,pad):
    """Returns dx,dy"""
    theta= np.arctan(abs(slope))
    return pad*np.sin(theta), pad*np.cos(theta)

def y1_line(rz,pad=None):
    slope,yint= 1.15,-0.15
    if pad:
        dx,dy= get_xy_pad(slope,pad)
        return slope*(rz+dx) + yint + dy
    else:
        return slope*rz + yint

def y2_line(rz,pad=None):
    slope,yint= -1.2,1.6
    if pad:
        dx,dy= get_xy_pad(slope,pad)
        return slope*(rz-dx) + yint + dy
    else:
        return slope*rz + yint

def get_ELG_box(rz,gr, pad=None):
    """
    Args:
        rz: r-z
        gr: g-r
        pad: magnitudes of padding to expand TS box
    """
    x1,y1= rz,y1_line(rz)
    x2,y2= rz,y2_line(rz)
    x3,y3= np.array([0.3]*len(rz)),gr
    x4,y4= np.array([1.6]*len(rz)),gr
    if pad:
        dx,dy= get_xy_pad(1.15,pad)
        x1,y1= x1-dx,y1+dy
        dx,dy= get_xy_pad(-1.2,pad)
        x2,y2= x2+dx,y2+dy
        x3 -= pad
        x4 += pad
    return dict(x1=x1, y1=y1,
                x2=x2, y2=y2,
                x3=x3, y3=y3,
                x4=x4, y4=y4)
#####

class TractorStamps(SimStamps):
    def __init__(self,ls_dir=None,outdir=None,
                 savedir=None, jpeg=False):
        """Same as SimStamps but for tractor catalogues

        Args:
            savedir: required for tractor not sims b/c cannot write to dr5 dir
        """
        super(TractorStamps,self).__init__(ls_dir=ls_dir,
                                           outdir=outdir,
                                           savedir=savedir,
                                           jpeg=jpeg)

    def set_paths_to_data(self,brick):
        """lists of catalogues filenames and coadd dirs"""
        self.cat_fns= [os.path.join(self.outdir,'tractor',
                                    brick[:3],
                                    'tractor-%s.fits' % brick)]
        self.coadd_dirs= [os.path.join(self.outdir,'coadd',
                                       brick[:3],brick)]
        if ((not os.path.exists(self.cat_fns[0])) |
            (not os.path.exists(self.coadd_dirs[0]))):
            raise OSError('does not exist: %s OR %s' % \
                    (self.cat_fns[0],self.coadd_dirs[0]))

    def set_xyid(self,zoom=None):
        x,y=self.cat.bx,self.cat.by
        if zoom:
            x -= zoom[0]
            y -= zoom[2]
        self.cat.set('x',x.astype(int))
        self.cat.set('y',y.astype(int))
        self.cat.set('id',self.cat.objid)

    def get_mag(self,band='g'):
        return flux2mag(self.cat.get('flux_'+band)/self.cat.get('mw_transmission_'+band))


    def apply_cuts(self):
        # Need extinction correction mag and colors
        d= {}
        for b in 'grz':
            d[b]= flux2mag(self.cat.get('flux_'+b)/self.cat.get('mw_transmission_'+b))
        df= pd.DataFrame(d)
        df['g-r']= df['g'] - df['r']
        df['r-z']= df['r'] - df['z']

        hasGRZ= ((self.cat.brick_primary) &
                 (self.cat.nobs_g >= 1) &
                 (self.cat.nobs_r >= 1) &
                 (self.cat.nobs_z >= 1))
        noArtifacts= ((self.cat.allmask_g == 0) &
                      (self.cat.allmask_r == 0) &
                      (self.cat.allmask_z == 0))

        keep= ((self.sim_sampling_cut(df)) &
               (self.isFaint_cut(df)) &
               #(noArtifacts) &
               (hasGRZ))

        len_bef= len(self.cat)
        self.cat.cut(keep)
        print('After cut, have %d/%d' % (len(self.cat),len_bef))


    def sim_sampling_cut(self,df):
        """same cut applied to simulated sources

        Args:
            df: pd.DataFrame have tractor cat extinction corrected grz mags
        """
        # TS box w/0.5 mag padding
        inBox= ((df['g-r'] <= y1_line(df['r-z'],pad=0.5)) &
                (df['g-r'] <= y2_line(df['r-z'],pad=0.5)) &
                (df['r-z'] >= 0.3 - 0.5) &
                (df['r-z'] <= 1.6 + 0.5))

        # Effective rhalf and drop comp, dev
        fwhm_or_rhalf= np.zeros(len(self.cat))-1 # arcsec
        isPSF= np.char.strip(self.cat.type) == 'PSF'
        isEXP= pd.Series(np.char.strip(self.cat.type)).isin(['EXP','REX'])
        isDEV= np.char.strip(self.cat.type) == 'DEV'
        isCOMP= np.char.strip(self.cat.type) == 'COMP'
        # rhalf ~ fwhm/2
        fwhm_or_rhalf[isPSF]= np.mean(np.array([self.cat[isPSF].psfsize_g,
                                                self.cat[isPSF].psfsize_r,
                                                self.cat[isPSF].psfsize_z]),axis=0)/2
        fwhm_or_rhalf[isEXP]= self.cat[isEXP].shapeexp_r
        fwhm_or_rhalf[isDEV]= self.cat[isDEV].shapedev_r

        grz_gt0= ((self.cat.flux_g > 0) &
                  (self.cat.flux_r > 0) &
                  (self.cat.flux_z > 0) &
                  (self.cat.flux_ivar_g > 0) &
                  (self.cat.flux_ivar_r > 0) &
                  (self.cat.flux_ivar_z > 0))

        keep= ((grz_gt0) &
               (isCOMP == False) &
               (isDEV == False) &
               (fwhm_or_rhalf < 5))

        #Last cut
        rhalf_lim= (0.262/2,2.) # Camera, Data
        g,r,z= tuple(np.array([24.0,23.4,22.5])+0.5)
        bad= ((fwhm_or_rhalf < rhalf_lim[0]) |
              (fwhm_or_rhalf > rhalf_lim[1]) |
              (df['z'] >= z) | #beyond mag limit
              (df['r'] >= r) |
              (df['g'] >= g))

        return (inBox) & (keep) & (bad == False)



    def isFaint_cut(self,df):
        """There are only faint sources in the deep2 matched sample,
        but in the tractor catalogus have a bright population presumably
        stars. Remove these

        Args:
            df: pd.DataFrame have tractor cat extinction corrected grz mags
        """
        # "elg_sample_5dim_10k.fits"
        min_mag= {'r': 20.4659,
                  'z': 19.4391,
                  'g': 20.6766}
        keep= ((df['g'] >= min_mag['g']) &
               (df['r'] >= min_mag['r']) &
               (df['z'] >= min_mag['z']))
        return keep

class UserDefinedStamps(SimStamps):
    def __init__(self,ls_dir=None,outdir=None,
                 savedir=None, jpeg=False):
        """Same as SimStamps but for tractor catalogues

        Args:
            savedir: required for tractor not sims b/c cannot write to dr5 dir
        """
        super().__init__(ls_dir=ls_dir,
                         outdir=outdir,
                         savedir=savedir,
                         jpeg=jpeg)

    def set_paths_to_data(self,brick):
        """lists of catalogues filenames and coadd dirs"""
        self.cat_fns= [os.path.join(self.savedir,'%s.fits' % brick)]
        self.coadd_dirs= [os.path.join(self.outdir,'coadd',
                                       brick[:3],brick)]
        if ((not os.path.exists(self.cat_fns[0])) |
            (not os.path.exists(self.coadd_dirs[0]))):
            raise OSError('does not exist: %s OR %s' % \
                    (self.cat_fns[0],self.coadd_dirs[0]))


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
             jpeg=False,
             bricks=[]):
    """

    Args:
        nproc: > 1 for mpi4py
        which: one of ['tractor','sim','userDefined']
        outdir: path to coadd,tractor dirs
        ls_dir: not needed if legacy_survey_dir env var already set
        savedir: where to write the hdf5 files, outdir if None
        jpeg: extract .jpg instead of .fits
        bricks: list bricks to make hdf5 cutouts from
    """
    assert(which in ['tractor','sim','userDefined'])
    if nproc > 1:
        from mpi4py.MPI import COMM_WORLD as comm
        bricks= np.array_split(bricks, comm.size)[comm.rank]

    d= dict(outdir=outdir,ls_dir=ls_dir,
            savedir=savedir, jpeg=jpeg)
    kwargs={}
    if which == 'sim':
        Obj= SimStamps(**d)
    elif which == 'tractor':
        Obj= TractorStamps(**d)
    elif which == 'userDefined':
        Obj= UserDefinedStamps(**d)
        kwargs.update(stampSize=64,
                      applyCuts=False)

    for brick in bricks:
        Obj.run(brick)


if __name__ == '__main__':
    #testcase_main()
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--which', type=str, choices=['tractor','sim','userDefined'], required=True, help='whether to make training hdf5 cutouts of real (DR5) or simulated (elg_dr5_coadds) outputs')
    parser.add_argument('--nproc', type=int, default=1, help='set to > 1 to run mpi4py')
    parser.add_argument('--bricks_fn', type=str, default=None, help='specify a fn listing bricks to run, or a single default brick will be ran')
    parser.add_argument('--savedir', type=str, default=None, help='specify a fn listing bricks to run, or a single default brick will be ran')
    parser.add_argument('--jpeg', action='store_true', default=False, help='put jpeg images in hdf5 file instead of fits coadss')
    args = parser.parse_args()

    # Data paths
    d= dict(savedir=args.savedir,
            jpeg=args.jpeg)
    if os.environ['HOME'] == '/home/kaylan':
        # ubuntu
        d.update(ls_dir=os.path.join(os.environ['HOME'],
                                 'mydata/legacysurveydir'))
        if args.which == 'sim':
            d.update(outdir=os.path.join(os.environ['HOME'],
                                 'mydata/elg_dr5_coadds'))
        elif args.which == 'tractor':
            d.update(outdir=os.path.join(os.environ['HOME'],
                                 'mydata/dr5_cutouts'))
            if args.savedir is None:
                d.update(savedir=os.path.join(os.environ['HOME'],
                                              'mydata/dr5_hdf5'))
    else:
        # nersc
        if args.which == 'sim':
            d.update(outdir=os.path.join(os.environ['CSCRATCH'],
                                 'obiwan_out/elg_dr5_coadds'))
        elif args.which == 'tractor':
            d.update(outdir='/global/project/projectdirs/cosmo/data/legacysurvey/dr5')
            if not args.savedir:
                d.update(savedir=os.path.join(os.environ['CSCRATCH'],
                                              'obiwan_out/dr5_hdf5'))
        elif args.which == 'userDefined':
            d.update(outdir='/global/project/projectdirs/cosmo/data/legacysurvey/dr5')
            if not args.savedir:
                d.update(savedir=os.path.join(os.environ['CSCRATCH'],
                                              'obiwan_out/dr5_hdf5'))

    # Bricks to run
    if not args.bricks_fn:
        bricks= ['1211p060'] #['1126p220']
    else:
        bricks= np.loadtxt(args.bricks_fn,dtype=str)

    mpi_main(nproc=args.nproc,which=args.which,bricks=bricks,
             **d)
