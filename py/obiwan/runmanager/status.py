"""
for a given brick, prints whether each rs* obiwan job finished or not
"""
import qdo
import os
import numpy as np
from glob import glob
from collections import defaultdict


QDO_RESULT= ['running', 'succeeded', 'failed']

def get_brickdir(outdir,obj,brick):
  return os.path.join(outdir,obj,brick[:3],brick)

def get_logdir(outdir,obj,brick,rowstart):
  return os.path.join( get_brickdir(outdir,obj,brick),
                       'rs%s' % rowstart)

def get_logfile(outdir,obj,brick,rowstart):
  return os.path.join( get_logdir(outdir,obj,brick,rowstart),
                       'log.%s' % brick)

def rs_from_logfile(log):
    logdir= os.path.dirname(log)
    return logdir[logdir.rfind('/rs')+3:]

def writelist(lis,fn):
  if os.path.exists(fn):
    os.remove(fn)
  with open(fn,'w') as foo:
    for li in lis:
      foo.write('%s\n' % li)
  print('Wrote %s' % fn)
  if len(lis) == 0:
    print('Warning: %s is empty list' % fn) 


class QdoList(object):
  def __init__(self,outdir,obj='elg',que_name='obiwan'):
    self.outdir= outdir
    self.obj= obj
    self.que_name= que_name
    print(self)

  def __str__(self):
    text= type(self).__name__ + ':' +'\n'
    for attr in ['outdir','obj','que_name']:
      text += '%s= %s\n' % (attr, getattr(self,attr) )
    return text

  def get_tasks_and_logs(self):
    """get tasks and logs for the three types of qdo status
    Running, Succeeded, Failed"""
    # Logs for all Failed tasks
    tasks={}
    logs={}
    #err= defaultdict(lambda: [])
    print('qdo Que: %s' % self.que_name)
    q = qdo.connect(self.que_name)
    for res in QDO_RESULT:
      # List of "brick rs" for each QDO_RESULT  
      tasks[res] = [a.task 
                    for a in q.tasks(state= getattr(qdo.Task, res.upper()))]
      # Corresponding log files  
      logs[res]= []  
      for task in tasks[res]:
        #brick,rs= np.array(task.split(' ')).astype(str) #Avoid unicode
        brick,rs= task.split(' ')
        logs[res].append( get_logfile(self.outdir,self.obj,brick,rs) )
    return tasks,logs


class RunStatus(object):
  """Tallys which QDO_RESULTS actually finished, what errors occured, etc.
  
  Args:
    tasks: dict, each key is list of qdo tasks
    logs: dict, each key is list of log files for each task
  """
    
  def __init__(self,tasks,logs):
    self.tasks= tasks
    self.logs= logs

  def tally(self):
    tally= defaultdict(list)
    if res == 'succeeded':
      for log in self.logs:
        with open(log,'r') as foo:
          text= foo.read()
        if "decals_sim:All done!" in text:
          tally[res].append( 1 )
        else:
          tally[res].append( 0 )

  def get_failed_errors(self):
      # Sort by type of error 
      if "EOFError: Ran out of input" in text:
          err['pickle'].append('%s %s' (brick,rs))
      elif "IndexError" in text:
          err['index'].append('%s %s' (brick,rs))
      elif "AssertionError" in text:
          err['assert'].append('%s %s' (brick,rs))

if __name__ == '__main__':
  from argparse import ArgumentParser
  parser = ArgumentParser()
  parser.add_argument('--outdir',default='/global/cscratch1/sd/kaylanb/obiwan_out/123/1238p245',help='',required=False)
  parser.add_argument('--obj',default='elg',help='',required=False)
  parser.add_argument('--brick',default='1238p245',help='',required=False)
  args = parser.parse_args()
  print(args)

  Q= QdoList(args.outdir,args.obj,que_name='obiwan_9deg')
  tasks,logs= Q.get_tasks_and_logs()
  
  # Write log fns so can inspect
  for res in logs.keys():
    writelist(logs[res],"%s_logfns.txt" % res)

  R= RunStatus(tasks,logs)
  #R.get_failed_errors()
  raise ValueError('done')

  cmd= os.path.join( get_brickdir(args.outdir,args.obj,args.brick) + '/*/' + 'log.*')
  logs= glob( cmd)
  assert(len(logs) > 0)

  status= defaultdict(dict)
  for log in logs:
    rs= rs_from_logfile(log)
    with open(log,'r') as foo:
      text= foo.read()
    if "INFO:decals_sim:All done!" in text:
      status[args.brick][rs]='done'
    else:
      status[args.brick][rs]='not'
  # Print results
  for brick in status.keys():
    for rs in status[brick].keys():
      print('%s %s: %s' % (brick,rs,status[brick][rs]))

  raise ValueError()



