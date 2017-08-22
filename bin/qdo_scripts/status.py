"""
Inspired by Martin Landriau LBNL
"""
import qdo
import os
from glob import glob
from collections import defaultdict

def get_log(outdir,obj,brick,rowstart):
  return os.path.join(outdir,obj,brick[:3],brick,'rs%d' % rowstart,
                      'log.%s' % brick)

def writelist(lis,fn):
  is os.path.exists(fn):
    os.remove(fn)
  with open(fn,'w') as foo:
    for li in lis:
      foo.write('%s\n' % li)
  print('Wrote %s' % fn)

outdir = os.path.join(os.environ['CSCRATCH'],
                      'obiwan_out/')
obj='elg'

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

