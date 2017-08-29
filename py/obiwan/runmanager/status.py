"""
for a given brick, prints whether each rs* obiwan job finished or not
"""
import qdo
import os
import numpy as np
from glob import glob
import re
from collections import defaultdict


QDO_RESULT= ['running', 'succeeded', 'failed']

def get_brickdir(outdir,obj,brick):
  return os.path.join(outdir,obj,brick[:3],brick)

def get_logdir(outdir,obj,brick,rowstart,doSkipid='no'):
  if doSkipid == "no":
    suffix= 'rs%s' % rowstart
  elif doSkipid == 'yes':
    suffix= 'skip_rs%s' % rowstart
  return os.path.join( get_brickdir(outdir,obj,brick),
                       '%s' % suffix)

def get_logfile(outdir,obj,brick,rowstart, doSkipid='no'):
  return os.path.join( get_logdir(outdir,obj,brick,rowstart, doSkipid=doSkipid),
                       'log.%s' % brick)

def get_slurm_files(outdir):
  return glob( outdir + '/slurm-*.out')

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
  """Queries the qdo db and maps log files to tasks and task status
  
  Args:
    outdir: obiwan outdir, the slurm*.out files are there
    obj: ...
    que_name: ie. qdo create que_name
  """

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

  def get_tasks_logs_slurms(self):
    """get tasks, logs, slurms for the three types of qdo status
    Running, Succeeded, Failed"""
    # Logs for all Failed tasks
    tasks={}
    ids={}
    logs= defaultdict(list)
    slurms= defaultdict(list)
    slurm_fns= get_slurm_files(self.outdir)
    assert(len(slurm_fns) > 0)
    print('slurm_fns=',slurm_fns)
    #err= defaultdict(lambda: [])
    print('qdo Que: %s' % self.que_name)
    q = qdo.connect(self.que_name)
    for res in QDO_RESULT:
      # List of "brick rs" for each QDO_RESULT  
      tasks[res] = [a.task 
                    for a in q.tasks(state= getattr(qdo.Task, res.upper()))]
      ids[res] = [a.id 
                    for a in q.tasks(state= getattr(qdo.Task, res.upper()))]
      # Corresponding log, slurm files  
      for task in tasks[res]:
        # Logs
        brick,rs,doSkipid = task.split(' ')
        logfn= get_logfile(self.outdir,self.obj,brick,rs, doSkipid=doSkipid)
        logs[res].append( logfn )
        # Slurms
        found= False
        for slurm_fn in slurm_fns:
          with open(slurm_fn,'r') as foo:
            text= foo.read()
          if logfn in text:
            found=True
            slurms[res].append( slurm_fn )
        if not found: 
            print('didnt find %s in slurms: ' % logfn,slurm_fn)
    return tasks,ids,logs,slurms

  def rerun_tasks(self,task_ids, debug=True):
    """set qdo tasks state to Pending for these task_ids
    
    Args:
      debug: False to actually reset the qdo tasks state AND to delete
      all output files for that task
    """
    q = qdo.connect(self.que_name)
    if not debug:
      print('resetting state for %d tasks' % len(task_ids))
    for task_id in task_ids:
      try:
        task_obj= q.tasks(id= task_id)
        brick,rs,doSkipid = task_obj.task.split(' ')
        logdir= get_logdir(self.outdir,self.obj,brick,rs, doSkipid=doSkipid)
        if debug:
          print('would remove dir %s' % logdir, 'for task obj',task_obj)
        else:
          task_obj.set_state(qdo.Task.PENDING)
          #os.removedirs(logdir)
      except ValueError:
        print('cant find task_id=%d' % task_id)

class RunStatus(object):
  """Tallys which QDO_RESULTS actually finished, what errors occured, etc.
  
  Args:
    tasks: dict, each key is list of qdo tasks
    logs: dict, each key is list of log files for each task

  Defaults:
    regex_errs: list of regular expressions matching possible log file errors
  """
    
  def __init__(self,tasks,logs):
    self.tasks= tasks
    self.logs= logs
    self.regex_errs= [
        'ValueError: starting row=[0-9]* exceeds number of artificial sources, quit',
        'No randoms in brick [0-9]*[pm][0-9]*, e.g. found nothing with db query:',
        'WARING: found nothing with:',
        'File "obiwan/kenobi.py", line 112, in get_skip_ids',
        'no skippedids.fits files exist for this brick'
        ]

  def get_tally(self):
    tally= defaultdict(list)
    for res in ['succeeded','failed']:
      if res == 'succeeded':
        for log in self.logs[res]:
          with open(log,'r') as foo:
            text= foo.read()
          if "decals_sim:All done!" in text:
            tally[res].append( 1 )
          else:
            tally[res].append( 0 )
      elif res == 'failed':
        for log in self.logs[res]:
          with open(log,'r') as foo:
            text= foo.read()
          found_err= False
          for regex in self.regex_errs:
            foundIt= re.search(regex, text)
            if foundIt:
              tally[res].append(regex)
              found_err=True
              break
          if not found_err:
            tally[res].append('Other')
    # numpy array, not list, works with np.where()
    for res in tally.keys():
      tally[res]= np.array(tally[res])
    return tally

  def print_tally(self,tally):
    for res in self.tasks.keys():
      print('--- Tally %s ---' % res)
      if res == 'succeeded':
         print('%d/%d = done' % (len(tally[res]), np.sum(tally[res])))
      elif res == 'failed':
        for regex in self.regex_errs + ['Other']:
          print('%d/%d = %s' % (
                   np.where(tally[res] == regex)[0].size, len(tally[res]), regex))
      elif res == 'running':
         print('%d/%d : need rerun' % (len(tally[res]),len(tally[res])))
  
  def get_logs_for_failed(self,regex='Other'):
    """Returns log and slurm filenames for failed tasks labeled as regex"""
    return self.logs[ tally['failed'] == regex ]



if __name__ == '__main__':
  from argparse import ArgumentParser
  parser = ArgumentParser()
  parser.add_argument('--qdo_quename',default='obiwan_9deg',help='',required=False)
  parser.add_argument('--outdir',default='/global/cscratch1/sd/kaylanb/obiwan_out/123/1238p245',help='',required=False)
  parser.add_argument('--obj',default='elg',help='',required=False)
  args = parser.parse_args()
  print(args)

  Q= QdoList(args.outdir,args.obj,que_name=args.qdo_quename)
  tasks,ids,logs,slurms= Q.get_tasks_logs_slurms()
  
  # Write log fns so can inspect
  for res in logs.keys():
    writelist(logs[res],"%s_%s_logfns.txt" % (args.qdo_quename,res))
    writelist(slurms[res],"%s_%s_slurmfns.txt" % (args.qdo_quename,res))

  R= RunStatus(tasks,logs)
  tally= R.get_tally()
  R.print_tally(tally)

  #err_logs= R.get_logs_for_failed(regex='Other')
  err_key= 'Other'
  err_logs= np.array(logs['failed'])[ tally['failed'] == err_key ]
  writelist(err_logs,"%s_%s_logsfailed.txt" % (args.qdo_quename,err_key))


  #Q.rerun_tasks(ids['running'], debug=False)
  raise ValueError('done')


