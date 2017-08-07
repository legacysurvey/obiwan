import fitsio
import os

def reduce_image(imgfn,newfn,extnames=['N4'],overwrite=False):
    h=fitsio.FITS(imgfn)
    if os.path.exists(newfn):
        if overwrite:
            os.remove(newfn)
        else:
            raise OSError('newfn=%s exists' % newfn)
    n=fitsio.FITS(newfn,'rw')
    n.write(None,header=h[0].read_header())
    for extnum in range(1,len(h)):
        if h[extnum].get_extname() in extnames:
            # legacypipe needs extnum and extname to match
            # so create intermediate empty hdus
            for new_extnum in range(1,extnum):
                n.write(np.array([-1]))
            n.write(h[extnum].read(),extname=h[extnum].get_extname())
    n.close()
    print('Wrote %s' % newfn)

if __name__ == "__main__":
    from astrometry.util.fits import fits_table, merge_tables
    import numpy as np
    import glob
    dr= '/project/projectdirs/desi/users/burleigh/obiwan/'
    a=fits_table(os.path.join(dr,'obiwan/elg/123/1238p245/rs0/coadd/legacysurvey-1238p245-ccds.fits'))
    for fn,ccdname in zip(np.char.strip(a.image_filename),np.char.strip(a.ccdname)):
        imgfn= os.path.join(dr,'legacy_survey_dir/images/',fn).replace('/decam/','/decam/DECam_CP/')
        newfn= imgfn.replace('/images/','/1ccd_images/')
        print(' imgfn=%s\n newfn=%s' % (imgfn,newfn))
        try:
            os.makedirs( os.path.dirname(newfn))
        except OSError:
            pass
        reduce_image(imgfn,newfn, extnames=[ccdname],overwrite=True)
        
    
