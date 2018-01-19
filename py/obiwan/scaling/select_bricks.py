"""
Using obiwan outputs, make a table of 'official' randoms per-brick
- uniform randoms: random ra,dec + geometry cut
- obiwan randoms: uniform randoms + recovered by tractor
"""

import numpy as np
import os
from glob import glob
import pandas as pd

try: 
    from astrometry.util.fits import fits_table, merge_tables
except ImportError:
    pass

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--data_dir', type=str, required=True, 
                        help='path to obiwan/, tractor/ dirs') 
    parser.add_argument('--bricks_with_rs3000_tractor', type=str, required=True, 
                        help='file listing bricks that have tractor cat in the rs3000 dir') 
    parser.add_argument('--outfn', type=str, required=False, 
                        default='brick_list.txt') 
    args = parser.parse_args()

    bricks= np.loadtxt(args.bricks_with_rs3000_tractor,dtype=str)
    np.random.seed(7)
    ind= np.random.randint(0,len(bricks),size=1000)
    bricks= np.array(list(set(bricks[ind])))

    keep=[]
    for brick in bricks:
        fn= os.path.join(args.data_dir,'coadd/%s/%s/rs3000' % (brick[:3],brick),
                         'legacysurvey-%s-ccds.fits' % brick)
        ccds= fits_table(fn)
        if len(set(ccds.filter)) == 3:
            print('found brick: %s' % brick)
            keep.append(brick)
        if len(keep) >= 10:
            break
    with open(args.outfn,'w') as foo:
        for kp in keep:
            foo.write('%s\n' % kp)
    print('Wrote %s' % args.outfn)


