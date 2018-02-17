"""
Images used in DR3 are now in different location on project/staging. This script
finds these cases and updated the survey-ccds files with the new imagefilename
"""

import numpy as np
import os
from glob import glob
import subprocess
from astrometry.util.fits import fits_table, merge_tables


def bashOutput(cmd):
    '''cmd is what would type on the command line, e.g. find . -name "hello.txt"'''
    proc = subprocess.Popen(cmd, shell=True,
                            stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    stdout, stderr= proc.communicate()
    return stdout,stderr
   
def mpi_main(nproc=1,imageFns=[]):
    """

    Args:
        nproc: > 1 for mpi4py
        bricks: list of bricks
    """
    if nproc > 1:
        from mpi4py.MPI import COMM_WORLD as comm
        imageFns= np.array_split(imageFns, comm.size)[comm.rank]
    old2new=[]
    notExist=[]
   
    dirs=dict(proj='/project/projectdirs/cosmo/staging/decam',
              proja='/global/projecta/projectdirs/cosmo/staging/decam')
    for cnt,imageFn in enumerate(imageFns):
        if cnt % 100 == 0:
            text= '%d/%d' % (cnt,len(imageFns))
            if nproc > 1:
                text= 'rank %d: ' % comm.rank + text
            print(text)
        if not os.path.exists(imageFn):
            found,err= bashOutput('find %s -name "%s"' % (dirs['proj'],os.path.basename(imageFn)))
            if len(found) == 0:
                found,err= bashOutput('find %s -name "%s"' % (dirs['proja'],os.path.basename(imageFn)))
            if len(found) == 0:
                #print("Does not exist anywhere: %s" % imageFn)
                notExist.append(imageFn)
            else:
                #print("--%s-- exists here --%s--" % (imageFn,found.decode().strip()))
                old2new.append((imageFn,found.decode().strip()))
    # Gather
    if nproc > 1:
        old2new = comm.gather(old2new, root=0)
        notExist = comm.gather(notExist, root=0)
    #if comm.rank == 0:
    fn='old2new.txt'
    with open(fn,'w') as foo:
        for i in range(len(old2new)):
            print(old2new[i][0],old2new[i][1])
            foo.write('%s --> %s\n' % (old2new[i][0],old2new[i][1]))
    print('Wrote %s' % fn)
    fn='notExist.txt'
    with open(fn,'w') as foo:
        for i in range(len(notExist)):
            foo.write('%s\n' % notExist[i])
    print('Wrote %s' % fn)
        

def imagefns_from_ccds_table(table_fn):
    T=fits_table(table_fn)
    fns= np.char.strip(T.image_filename)
    return np.array(list(set(fns)))


if __name__ == '__main__':
    #testcase_main()
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--nproc', type=int, default=1, help='set to > 1 to run mpi4py') 
    parser.add_argument('--image_list', type=str, default=None, 
                        help='specify a fn listing the image filenames') 
    parser.add_argument('--ccds_table', type=str, default=None, 
                        help='use a ccd table insteady, grabbing the imagefilename column') 
    args = parser.parse_args()
    
    if args.image_list:
        imageFns= np.loadtxt(args.image_list,dtype=str)
    elif args.ccds_table:
        imageFns= imagefns_from_ccds_table(args.ccds_table)

    mpi_main(nproc=args.nproc,imageFns=imageFns)

