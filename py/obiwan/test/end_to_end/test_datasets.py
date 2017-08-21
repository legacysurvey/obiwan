# tests datasets DR3, DR5, DR3_eBOSS

from __future__ import print_function
import os

from obiwan.kenobi import main,get_parser

DATASETS= ['DR3','DR5','DR3_eBOSS']

def run_dataset(dataset):
    """run a single dataset's test problem

    Args:
      dataset: string, 'DR3', 'DR5', etc
    """
    assert(dataset in DATASETS)
    print('testing dataset: %s' % dataset)
    obj='elg'
    brick='1238p245'
    os.environ["LEGACY_SURVEY_DIR"]= os.path.join(os.path.dirname(__file__), 
                                                  'legacypipedir_%s_dataset_%s' % (brick,dataset))
    outdir = os.path.join(os.path.dirname(__file__),
                          'outdir_%s_%s' % (brick,dataset))
    randoms_from_fits= os.path.join(os.path.dirname(__file__), 
                                    'randoms/%s' % obj,'randoms_%s_10.fits' % brick)

    cmd_line=['--dataset', dataset, '-b', brick, '-n', '2', 
              '-o', 'elg', '--outdir', outdir, '--add_sim_noise',
              '--randoms_from_fits', randoms_from_fits]
    parser= get_parser()
    args = parser.parse_args(args=cmd_line)
    
    main(args=args)
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

def test_datasets():
    run_dataset('DR3')
    run_dataset('DR5')
    assert(True)

if __name__ == "__main__":
    run_dataset('DR5')
    

    
