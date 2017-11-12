from __future__ import print_function
import os
import sys
import numpy as np
from legacypipe.runbrick import main
#from legacyanalysis.decals_sim import main as sim_main
from astrometry.util.fits import fits_table

os.environ["DUST_DIR"]= os.path.join(os.environ['HOME'],'myrepo/dust')

# RexGalaxy, DECam data
surveydir = os.path.join(os.path.dirname(__file__), 'testcase6')
outdir = os.path.join(os.path.dirname(__file__), 'testcase6-out')
print('surveydir=%s' % surveydir)
main(args=['--brick', '1102p240', '--zoom', '500', '600', '650', '750',
           '--force-all', '--no-write', '--no-wise',
           #'--rex', #'--plots',
           '--survey-dir', surveydir,
           '--outdir', outdir])
fn = os.path.join(outdir, 'tractor', '110', 'tractor-1102p240.fits')
assert(os.path.exists(fn))
T = fits_table(fn)
assert(len(T) == 2)
print('Types:', T.type)
assert(T.type[0] == 'REX ')
cmd = ('(cd %s && sha256sum -c %s)' %
       (outdir, os.path.join('tractor', '110', 'brick-1102p240.sha256sum')))
print(cmd)
rtn = os.system(cmd)
assert(rtn == 0)
    

    
