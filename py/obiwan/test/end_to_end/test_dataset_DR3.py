from __future__ import print_function
import os

from obiwan.kenobi import main

#os.environ["DUST_DIR"]= os.path.join(os.environ['HOME'],'myrepo/dust')

dataset='DR3'
brick='1238p245'
os.environ["LEGACY_SURVEY_DIR"]= os.path.join(os.path.dirname(__file__), 
                                'legacypipedir_%s_dataset_%s' % (brick,dataset))
outdir = os.path.join(os.path.dirname(__file__), 'outdir_%s_%s' % (brick,dataset))

main(args=['--dataset', dataset, '-b', brick, '-n', '2', 
           '-o', 'elg', '--outdir', outdir, '--add_sim_noise'])
assert(True)
#fn = os.path.join(outdir, 'tractor', '110', 'tractor-1102p240.fits')
#assert(os.path.exists(fn))
#T = fits_table(fn)
#assert(len(T) == 2)
#print('Types:', T.type)
#assert(T.type[0] == 'REX ')
#cmd = ('(cd %s && sha256sum -c %s)' %
#       (outdir, os.path.join('tractor', '110', 'brick-1102p240.sha256sum')))
#print(cmd)
#rtn = os.system(cmd)
#assert(rtn == 0)
    

    
