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
   
def mpi_fns_that_exist(nproc=1,imageFns=[],suffix=''):
    """

    Args:
        nproc: > 1 for mpi4py
        bricks: list of bricks
        suffix: identifying name to append to output fns
    """
    if nproc > 1:
        from mpi4py.MPI import COMM_WORLD as comm
        imageFns= np.array_split(imageFns, comm.size)[comm.rank]
    else:
        class MyComm(object):
            def __init__(self):
                self.rank=0
        comm= MyComm()
    old2new=[]
    notExist=[]
   
    dirs=dict(proj='/project/projectdirs/cosmo/staging',
              proja='/global/projecta/projectdirs/cosmo/staging')
    for cnt,imageFn in enumerate(imageFns):
        if cnt % 100 == 0:
            text= 'rank %d: %d/%d' % (comm.rank,cnt,len(imageFns))
            print(text)
        onProj= os.path.join(dirs['proj'],imageFn)
        onProja= os.path.join(dirs['proja'],imageFn)
        if os.path.exists(onProj):
            pass
        elif os.path.exists(onProja): 
            pass
        else:
            print('Doesnt exist on proj or proja under orig name: %s' % imageFn)
            found,err= bashOutput('find %s -name "%s"' % (dirs['proj'],os.path.basename(imageFn)))
            if len(found) == 0:
                found,err= bashOutput('find %s -name "%s"' % (dirs['proja'],os.path.basename(imageFn)))
            if len(found) == 0:
                print("Doesnt exist anywhere: %s" % imageFn)
                notExist.append(imageFn)
            else:
                newfn= (found.decode().strip()
                             .replace(dirs['proj'],'')
                             .replace(dirs['proja'],''))
                old2new.append((imageFn,newfn))
        else:
            print('Exists: %s' % imageFn)
    # Gather
    if nproc > 1:
        old2new = comm.gather(old2new, root=0)
        notExist = comm.gather(notExist, root=0)
    #if comm.rank == 0:
    #if nproc > 1:
    #    fn= fn.replace('.txt','_rank%d.txt' % comm.rank)
    if comm.rank == 0: 
        fn='old2new%s.txt' % suffix
        with open(fn,'w') as foo:
            for i in range(len(old2new)):
                print(old2new[i][0],old2new[i][1])
                foo.write('%s --> %s\n' % (old2new[i][0],old2new[i][1]))
        print('Wrote %s' % fn)
        
        fn='notExist%s.txt' % suffix
        #if nproc > 1:
        #    fn= fn.replace('.txt','_rank%d.txt' % comm.rank)
        with open(fn,'w') as foo:
            for i in range(len(notExist)):
                foo.write('%s\n' % notExist[i])
        print('Wrote %s' % fn)
        

def imagefns_from_ccds_table(table_fn):
    T=fits_table(table_fn)
    fns= np.char.strip(T.image_filename)
    return np.array(list(set(fns)))


def update_surveyccds(ccds_fn,old2new_fn):
    ccds=fits_table(ccds_fn)
    saveName= (ccds_fn.replace('.fits.gz','-newfns.fits.gz')
                      .replace('.fits','-newfns.fits'))
    with open(old2new_fn,'r') as f:
        data= f.readlines()
    old=[row.split(" ")[0] for row in data]
    new=[row.split(" ")[-1].strip() for row in data]
    
    fns= np.char.strip(ccds.image_filename)
    for o,n in zip(old,new):
        fns[fns == o]= n
    ccds.set('image_filename',fns)
    ccds.writeto(saveName)
    print('Wrote %s' % saveName)




if __name__ == '__main__':
    #testcase_main()
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--runWhat', type=str, choices=['fns_that_exist','update_ccds'],required=True) 
    parser.add_argument('--nproc', type=int, default=1, help='set to > 1 to run mpi4py') 
    parser.add_argument('--image_list', type=str, default=None, 
                        help='specify a fn listing the image filenames') 
    parser.add_argument('--ccds_table', type=str, default=None, 
                        help='use a ccd table insteady, grabbing the imagefilename column') 
    parser.add_argument('--old2new_fn', type=str, default=None, 
                        help='output by mpi_fns_that_exist') 
    args = parser.parse_args()
    
    if args.runWhat == 'fns_that_exist':
        if args.image_list:
            imageFns= np.loadtxt(args.image_list,dtype=str)
            suffix= os.path.basename(args.image_list).replace('.txt','')
        elif args.ccds_table:
            imageFns= imagefns_from_ccds_table(args.ccds_table)
            suffix= os.path.basename(args.ccds_table).replace('.fits.gz','').replace('.fits','')
        mpi_fns_that_exist(nproc=args.nproc,imageFns=imageFns,suffix='_'+suffix)
    elif args.runWhat == 'update_ccds':
        update_surveyccds(args.ccds_table,args.old2new_fn)
