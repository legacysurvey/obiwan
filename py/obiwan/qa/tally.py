import os
import pandas as pd
from glob import glob

from obiwan.runmanager.status import QdoList,writelist
from obiwan.kenobi import dobash

def bricksFromLogs(loglist):
  """return set of bricks occuring in logfiles
  
  Args:
    loglist: list of filenames for log files
  """
  logs= pd.DataFrame({'logs':loglist})
  bricks= logs.loc[:,'logs'].str.split('log.').str[-1]
  return list(set(bricks))

def bricksInDR5(bricklist):
  """return bricks in bricklist that are in DR5"""
  dr5_out='/global/cscratch1/sd/desiproc/DR5_out'
  dr5_bricks=[] 
  for brick in bricklist:
    tractor_fn= os.path.join(dr5_out,'tractor',brick[:3],
                             'tractor-%s.fits' % brick)
    if os.path.exists(tractor_fn):
      dr5_bricks.append( brick )
  inBoth= set(bricklist).intersection(set(dr5_bricks))
  return list(inBoth)


if __name__ == '__main__':
  from argparse import ArgumentParser
  parser = ArgumentParser()
  parser.add_argument('--qdo_quename',default='obiwan_9deg',help='',required=False)
  parser.add_argument('--outdir',default='/global/cscratch1/sd/kaylanb/obiwan_out/123/1238p245',help='',required=False)
  parser.add_argument('--obj',default='elg',help='',required=False)
  args = parser.parse_args()
  print(args)
  
  Q= QdoList(args.outdir,args.obj,que_name=args.qdo_quename)
  result= 'succeeded'
  tasks,ids,logs_dict,slurms= Q.get_tasks_logs_slurms(one_result= result, getslurms=False)
  logs= logs_dict[result]
  writelist(logs,"%s_%s_logfns.txt" % (args.qdo_quename,result))
  
  bricks= bricksFromLogs(logs)
  bricks_both= bricksInDR5(bricks)
  if len(bricks_both) == 0:
    raise ValueError('none of bricks are in DR5!')
  writelist(bricks_both,"qa_bricksInBoth.txt")


  


