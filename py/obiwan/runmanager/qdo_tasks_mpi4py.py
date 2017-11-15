from __future__ import division, print_function
from argparse import ArgumentParser
import numpy as np

from obiwan.runmanager.qdo_tasks import TaskList

parser = ArgumentParser()
parser.add_argument('--ra1',type=float,action='store',help='bigbox',required=True)
parser.add_argument('--ra2',type=float,action='store',help='bigbox',required=True)
parser.add_argument('--dec1',type=float,action='store',help='bigbox',required=True)
parser.add_argument('--dec2',type=float,action='store',help='bigbox',required=True)
parser.add_argument('--nobj_per_run',type=float,action='store',help='bigbox',required=True)
parser.add_argument('--obj', type=str, choices=['star','elg'], default=None, required=True) 
parser.add_argument('--randoms_db',type=str,action='store',required=True)
parser.add_argument('--outdir',type=str,action='store',help='/global/cscratch1/sd/kaylanb/obiwan_out/elg_9deg2_ra175',required=True)
parser.add_argument('--do_more',type=str,action='store',default='no',required=False)
parser.add_argument('--minid',type=int,action='store',default=1,required=False)
args = parser.parse_args()

T= TaskList(ra1=args.ra1,ra2=args.ra2, dec1=args.dec1, dec2=args.dec2,
            nobj_per_run=args.nobj_per_run)
do_skipids='no'
T.get_bricks()

# split bricks
from mpi4py.MPI import COMM_WORLD as comm
bricks= np.array_split(T.bricks,comm.size)[comm.rank] 
tasks, not_in_bricks= T.get_tasklist(objtype=args.obj,randoms_db=args.randoms_db,
                                     do_more=args.do_more,minid=args.minid,
                                     outdir=args.outdir,
                                     bricks=bricks)
# pool tasks and not_in_bricks lists
all_tasks = comm.gather(tasks,root=0)
all_not_bricks = comm.gather(not_in_bricks,root=0)
T.writetasks(all_tasks,all_not_bricks,
             do_skipids,args.do_more,args.minid)


