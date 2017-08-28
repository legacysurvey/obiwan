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

DOWNLOAD_ROOT = "http://portal.nersc.gov/project/desi/users/kburleigh/obiwan/"
NERSC_ROOT = DOWNLOAD_ROOT.replace("http://portal.nersc.gov/project/",
                                   "/global/project/projectdirs/")\
                          .replace("/users/","/www/users/")

class TestData(object):    
    def fetch(self,outdir):
        name= 'testdata.tar.gz'
        self.outdir= os.path.join(outdir,name.replace('.tar.gz',''))
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
        remote_fn= DOWNLOAD_ROOT + name
        local_fn= os.path.join(outdir,name)
        if not os.path.exists(local_fn):
            print('Retrieving %s, extracting here %s' % (remote_fn,local_fn))
            urllib.request.urlretrieve(remote_fn, local_fn)
            tgz = tarfile.open(local_fn)
            tgz.extractall(path=outdir)
            tgz.close()
        
    def get_fn(self,name):
        if name == 'skipids':
            fn= 'skipids_table.fits'
        else:
            raise ValueError('name=%s not supported' % name)
        return os.path.join(self.outdir,fn)
    

    

    
