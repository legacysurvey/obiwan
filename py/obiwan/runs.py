"""
Grabs, configures, and otherwise sets up a production run on eBOSS DR3 data or DR5 data
"""

import os
from glob import glob
import numpy as np

from theValidator.catalogues import CatalogueFuncs

DOWNLOAD_ROOT = "http://portal.nersc.gov/project/desi/users/kburleigh/obiwan/"
NERSC_ROOT = DOWNLOAD_ROOT.replace("http://portal.nersc.gov/project/",
                                   "/global/project/projectdirs/")\
                          .replace("/users/","/www/users/")

def in_eboss(T):
    x= np.logical_and.reduce((T.ra > 317.,T.ra < 360., T.dec > -2, T.dec < 2))
    y= np.logical_and.reduce((T.ra > 0.,T.ra < 45., T.dec > -5.,T.dec < 5.))
    z= np.logical_and.reduce((T.ra > 126., T.ra < 169.,T.dec > 14.,T.dec < 29.))
    return np.logical_or.reduce((x,y,z)) 

if __name__ == "__main__":
  eboss_fns= glob(os.path.join(NERSC_ROOT,'configs/dr3eBOSS/additional_ccds',
                               "survey-ccds-*.fits.gz"))
  dr3_fns= glob(os.path.join(NERSC_ROOT,'configs/dr3eBOSS/dr3',
                             "survey-ccds-*.fits.gz"))


  dr3= CatalogueFuncs().stack(dr3_fns,textfile=False)
  eboss= CatalogueFuncs().stack(eboss_fns,textfile=False)

  dr3.set('pid', np.char.add(dr3.expnum.astype(str),
                             np.char.strip(dr3.image_filename)))
  eboss.set('pid', np.char.add(eboss.expnum.astype(str),
                               np.char.strip(eboss.image_filename))

  dr3.cut( in_eboss(dr3))
  eboss.cut( in_eboss(eboss))

  len( set(dr3.pid).difference(set(eboss.pid)) )
  len( set(eboss.pid).difference(set(dr3.pid)) )
  len( set(eboss.pid).intersection(set(dr3.pid)) )
