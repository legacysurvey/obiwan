"""
for a given brick, prints whether each rs* obiwan job finished or not
"""
import qdo
import os
import numpy as np
from glob import glob
import re
from collections import defaultdict

from obiwan.common import get_savedir,dobash

QDO_RESULT= ['running', 'succeeded', 'failed']


def get_logfile(outdir,obj,brick,rowstart, 
                do_skipids='no',do_more='no'):
  return os.path.join( get_savedir(outdir,obj,brick,rowstart, 
                                   do_skipids=do_skipids,do_more=do_more),
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

  def get_tasks_logs(self):
    """get tasks, logs, slurms for the three types of qdo status
    Running, Succeeded, Failed
    
    Args:
      one_result: only get logs slurms etc for one of succeeded, failed, running
      getslurms: set to False to skip finding slurm.out file
    """
    # Logs for all Failed tasks
    tasks={}
    ids={}
    logs= defaultdict(list)
    #err= defaultdict(lambda: [])
    print('qdo Que: %s' % self.que_name)
    q = qdo.connect(self.que_name)
    if one_result:
      assert(one_result in QDO_RESULT)
      all_results= [one_result]
    else:
      all_results= QDO_RESULT
    for res in all_results:
      # List of "brick rs" for each QDO_RESULT  
      tasks[res] = [a.task 
                    for a in q.tasks(state= getattr(qdo.Task, res.upper()))]
      ids[res] = [a.id 
                    for a in q.tasks(state= getattr(qdo.Task, res.upper()))]
      # Corresponding log, slurm files  
      for task in tasks[res]:
        # Logs
        brick,rs,do_skipids,do_more = task.split(' ')
        logfn= get_logfile(self.outdir,self.obj,brick,rs, 
                           do_skipids=do_skipids,do_more=do_more)
        logs[res].append( logfn )
    return tasks,ids,logs

  def rerun_tasks(self,task_ids, modify=False):
    """set qdo tasks state to Pending for these task_ids
    
    Args:
      modify: True to actually reset the qdo tasks state AND to delete
        all output files for that task
    """
    q = qdo.connect(self.que_name)
    for task_id in task_ids:
      try:
        task_obj= q.tasks(id= int(task_id))
        brick,rs,do_skipids,do_more= task_obj.task.split(' ')
        logdir= get_savedir(self.outdir,self.obj,brick,rs, 
                            do_skipids=do_skipids,do_more=do_more)
        rmcmd= "rm %s/*" % logdir
        if modify:
          task_obj.set_state(qdo.Task.PENDING)
          dobash(rmcmd)
        else:
          print('would remove id=%d, which corresponds to taks_obj=' % task_id,task_obj)
          print('by calling dobash(%s)' % rmcmd)
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
        'no skippedids.fits files exist for this brick',
        "KeyError: 'do_skipds'",
        'psycopg2.OperationalError: FATAL:  remaining connection slots are reserved',
        'psycopg2.OperationalError: FATAL:  sorry, too many clients already']
    self.regex_errs_extra= ['Other','log does not exist']

  def get_tally(self):
    tally= defaultdict(list)
    for res in ['succeeded','failed','running']:
      if res == 'succeeded':
        for log in self.logs[res]:
          with open(log,'r') as foo:
            text= foo.read()
          if "decals_sim:All done!" in text:
            tally[res].append( 1 )
          else:
            tally[res].append( 0 )
      elif res == 'running':
        for log in self.logs[res]:
          tally[res].append(1)
      elif res == 'failed':
        for log in self.logs[res]:
          if not os.path.exists(log):
            tally[res].append('log does not exist')
            continue
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
        for regex in self.regex_errs + self.regex_errs_extra:
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
  parser.add_argument('--outdir',default='/global/cscratch1/sd/kaylanb/obiwan_out/elg_9deg2_ra175',help='',required=False)
  parser.add_argument('--obj',default='elg',help='',required=False)
  parser.add_argument('--running_to_pending',action="store_true",default=False,help='set to reset all "running" jobs to "pending"')
  parser.add_argument('--failed_message_to_pending',action='store',default=None,help='set to message of failed tak and reset all failed tasks with that message to pending')
  parser.add_argument('--modify',action='store_true',default=False,help='set to actually reset the qdo tasks state AND to delete IFF running_to_pending or failed_message_to_pending are set')
  args = parser.parse_args()
  print(args)

  Q= QdoList(args.outdir,args.obj,que_name=args.qdo_quename)
  tasks,ids,logs= Q.get_tasks_logs()
  
  # Write log fns so can inspect
  for res in logs.keys():
    writelist(logs[res],"%s_%s_logfns.txt" % (args.qdo_quename,res))

  R= RunStatus(tasks,logs)
  tally= R.get_tally()
  R.print_tally(tally)

  # logs,tasks for each type of failure
  for err_key in R.regex_errs + R.regex_errs_extra:
    err_logs= np.array(logs['failed'])[ tally['failed'] == err_key ]
    err_tasks= np.array(tasks['failed'])[ tally['failed'] == err_key ]
    err_string= (err_key[:12] + err_key[-8:])
    for rem_str in [" ",":",",","*",'"']:
      err_string= err_string.replace(rem_str,"_")
    writelist(err_logs,"logs_%s_%s.txt" % (args.qdo_quename,err_string))
    writelist(err_tasks,"tasks_%s_%s.txt" % (args.qdo_quename,err_string))

  # Rerun tasks and delete those tasks' outputs
  if args.running_to_pending:
    if len(ids['running']) > 0:
      Q.rerun_tasks(ids['running'], modify=args.modify)
  if args.failed_message_to_pending:
    hasMessage= np.where(tally['failed'] == args.failed_message_to_pending)[0]
    if hasMessage.size > 0:
      theIds= np.array(ids['failed'])[hasMessage]
      Q.rerun_tasks(theIds, modify=args.modify)


  print('done')


