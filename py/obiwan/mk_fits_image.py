import fitsio
import os

def reduce_image(imgfn,newfn,extname='N4',overwrite=False):
    print(' imgfn=%s\n newfn=%s' % (imgfn,newfn))
    h=fitsio.FITS(imgfn)
    extname= extname.upper()
    extnum= h[extname].get_extnum()
    if os.path.exists(newfn):
        if overwrite:
            os.remove(newfn)
        else:
            raise OSError('newfn=%s exists' % newfn)
    n=fitsio.FITS(newfn,'rw')
    n.write(None,header=h[0].read_header())
    # legacypipe needs extnum and extname to match
    # so create intermediate empty hdus
    for junk in range(1,extnum):
        n.write(np.array([-1]))
    n.write(h[extnum].read(),header=h[extnum].read_header(),
            extname=h[extnum].get_extname())
    assert(n[extnum].get_extnum() == extnum)
    n.close()
    print('Wrote %s' % newfn)

def reduce_images(imgfn,newfn,**kwargs):
    reduce_image(imgfn,newfn,**kwargs)
    reduce_image(imgfn.replace('ooi','oow'),newfn.replace('ooi','oow'),**kwargs)
    reduce_image(imgfn.replace('ooi','ood'),newfn.replace('ooi','ood'),**kwargs)

if __name__ == "__main__":
    from astrometry.util.fits import fits_table, merge_tables
    import numpy as np
    import glob
    dr= '/project/projectdirs/desi/users/burleigh/obiwan/'
    a=fits_table(os.path.join(dr,'legacy_survey_dir/survey-ccds-decals-7ccds.fits.gz'))
    for fn,ccdname in zip(np.char.strip(a.image_filename),np.char.strip(a.ccdname)):
        imgfn= os.path.join(dr,'legacy_survey_dir/allccd_images/',fn).replace('/decam/','/decam/DECam_CP/')
        newfn= imgfn.replace('/allccd_images/','/1ccd_images/')
        try:
            os.makedirs( os.path.dirname(newfn))
        except OSError:
            pass
        kwargs= dict(extname=ccdname, overwrite=True)
        reduce_images(imgfn,newfn, **kwargs)
        
    
