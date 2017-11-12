# See LICENSE.rst for BSD 3-clause license info
# -*- coding: utf-8 -*-
"""
Fetches and loads test data that is too large for github
"""

import fitsio
from astrometry.util.fits import fits_table, merge_tables
import os
import sys
from glob import glob

from obiwan.fetch import fetch_targz

DOWNLOAD_DIR='http://portal.nersc.gov/project/desi/users/kburleigh/legacyzpts'

class GetData(object):  
    def fetch(self):
        name= 'testdata.tar.gz'
        fetch_targz(os.path.join(DOWNLOAD_DIR,                                             
                                'ccds_decam.tar.gz'),                                
                    self.outdir())

    def load(self):
        fns= glob(os.path.join(self.outdir(),'ccds_decam',
                  'small_c4d*.fits.fz'))
        hdu= fitsio.FITS(fns[0])
        return hdu[1].read()

    def outdir(self):
        return os.path.join(os.path.dirname(__file__),
                            'testccds')

    
if __name__ == '__main__':
    d= GetData()
    d.fetch()
    img= d.load()

    

    
