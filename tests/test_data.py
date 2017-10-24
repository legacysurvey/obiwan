# See LICENSE.rst for BSD 3-clause license info
# -*- coding: utf-8 -*-
"""
Fetches and loads test data that is too large for github
"""

from astrometry.util.fits import fits_table, merge_tables
import os
import sys
from glob import glob
from six.moves import urllib
import tarfile

from obiwan.fetch import fetch_targz


DOWNLOAD_ROOT = "http://portal.nersc.gov/project/desi/users/kburleigh/obiwan/"
NERSC_ROOT = DOWNLOAD_ROOT.replace("http://portal.nersc.gov/project/",
                                   "/global/project/projectdirs/")\
                          .replace("/users/","/www/users/")

class TestData(object):  
    def fetch(self,outdir):
        name= 'testdata.tar.gz'
        fetch_targz(os.path.join(DOWNLOAD_ROOT,name),
                    outdir)
           
    def get_fn(self,name,outdir):
        if name == 'skipids':
            fn= 'skipids_table.fits'
        else:
            raise ValueError('name=%s not supported' % name)
        return os.path.join(outdir,fn)
    

    

    
