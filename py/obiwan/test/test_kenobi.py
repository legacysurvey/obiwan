import os

from astrometry.util.fits import fits_table, merge_tables
from astrometry.libkd.spherematch import match_radec
from obiwan import kenobi
from obiwan.common import get_savedir
from obiwan.test.test_data import TestData

def test_output_dir():
    kwargs=dict(outdir= 'hello',
                obj='elg',
                brick='1238p245',
                rowstart=1,
                do_skipids='no',
                do_more='no')
    ans= 'hello/elg/123/1238p245/rs1'
    assert(kenobi.get_savedir(**kwargs) == ans)
    kwargs.update(do_skipids='yes')
    ans= 'hello/elg/123/1238p245/skip_rs1'
    assert(get_savedir(**kwargs) == ans)
    kwargs.update(do_more='yes')
    ans= 'hello/elg/123/1238p245/more_skip_rs1'
    assert(get_savedir(**kwargs) == ans)

def test_overlapping_radec():
    data= TestData()
    data.fetch(outdir= os.path.dirname(__file__))
    fn= data.get_fn('skipids')
    Samp= fits_table(fn)
    for iter in range(5):
        keep,rem= kenobi.flag_nearest_neighbors(Samp, radius_in_deg=5./3600)
        #print('iter %d: Samp size %d, number close sources %d' % (iter,len(Samp),len(rem)))
        Samp.cut(keep)
    # No more nearby sources in the sample
    assert(len(rem) == 0)                     
    
