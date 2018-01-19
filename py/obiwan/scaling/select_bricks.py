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

def ten_random_bricks(bricks,data_dir,outfn):
    np.random.seed(7)
    ind= np.random.randint(0,len(bricks),size=1000)
    bricks= np.array(list(set(bricks[ind])))

    keep=[]
    for brick in bricks:
        fn= os.path.join(data_dir,'coadd/%s/%s/rs3000' % (brick[:3],brick),
                         'legacysurvey-%s-ccds.fits' % brick)
        ccds= fits_table(fn)
        if len(set(ccds.filter)) == 3:
            print('found brick: %s' % brick)
            keep.append(brick)
        if len(keep) >= 10:
            break
    with open(outfn,'w') as foo:
        for kp in keep:
            foo.write('%s\n' % kp)
    print('Wrote %s' % outfn)

def qdo_tasks(bricks,nobj=500):
    outfn='qdo_tasks_nobj%d.txt' % nobj
    with open(outfn,'w') as foo:
        for brick in bricks:
            num=0
            while(num < 1500):
                foo.write('%s %d no no\n' % (brick,num))
                num+= nobj
                if num >= 1500:
                    break
    print('Wrote %s' % outfn)

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--data_dir', type=str, required=True, 
                        help='path to obiwan/, tractor/ dirs') 
    parser.add_argument('--bricks_with_rs3000_tractor', type=str, required=True, 
                        help='file listing bricks that have tractor cat in the rs3000 dir') 
    parser.add_argument('--ten_bricks_fn', type=str, default='ten_bricks.txt',required=False, 
                        help='fn of list of 10 bricks to make qdo tasks for') 
    args = parser.parse_args()


    bricks= np.loadtxt(args.bricks_with_rs3000_tractor,dtype=str)
    if not os.path.exists(args.ten_bricks_fn):
        ten_random_bricks(bricks,args.data_dir,args.ten_bricks_fn)
    ten_bricks= np.loadtxt(args.ten_bricks_fn,dtype=str)
    for n in [500,1000,1500]:
        qdo_tasks(ten_bricks,nobj=n)


