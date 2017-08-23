"""
for a given brick, prints whether each rs* obiwan job finished or not
"""
import qdo
import os
from glob import glob
from collections import defaultdict

def get_brickdir(outdir,obj,brick):
  return os.path.join(outdir,obj,brick[:3],brick)

def get_logdir(outdir,obj,brick,rowstart):
  return os.path.join( get_brickdir(outdir,obj,brick),
                       'rs%d' % rowstart)

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


if __name__ == '__main__':
  from argparse import ArgumentParser
  parser = ArgumentParser()
  parser.add_argument('--brick',default='1238p245',help='',required=False)
  parser.add_argument('--outdir',default='/global/cscratch1/sd/kaylanb/obiwan_out/123/1238p245',help='',required=False)
  parser.add_argument('--obj',default='elg',help='',required=False)
  args = parser.parse_args()

  print(args)
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
  # Logs for all Failed tasks
  tasks={}
  logs_fail=[]
  err= defaultdict(lambda: [])
  q = qdo.connect('obiwan')
  tasks['finish'] = q.tasks(state=qdo.Task.SUCCEEDED)
  tasks['fail'] = q.tasks(state=qdo.Task.FAILED)
  for i in range(len(tasks['fail'])):
      brick,rs= tasks['fail'][i]
      logfn=  get_log(outdir,obj,brick,rs)
      logs_fail.append(logfn)
      # Sort by type of error 
      with open(logfn,'r') as foo:
        text= foo.read()
      if "EOFError: Ran out of input" in text:
          err['pickle'].append('%s %s' (brick,rs))
      elif "IndexError" in text:
          err['index'].append('%s %s' (brick,rs))
      elif "AssertionError" in text:
          err['assert'].append('%s %s' (brick,rs))

  # Write results
  for key in err.keys():
    if err[key].size > 0:
      writelist(err[key],"%s.txt" % key)

  print(logs_fail)

