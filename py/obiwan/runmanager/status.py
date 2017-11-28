"""
for a given brick, prints whether each rs* obiwan job finished or not
"""
import os
import numpy as np
from glob import glob
import re
from collections import defaultdict

from obiwan.common import dobash,writelist,get_rsdir

import qdo

QDO_RESULT= ['running', 'succeeded', 'failed']


def get_interm_dir(outdir,brick,rowstart, 
                do_skipids='no',do_more='no'):
    """Returns paths like outdir/bri/brick/rs0"""
    rsdir= get_rsdir(rowstart, 
                     do_skipids,do_more)
    return os.path.join(outdir,brick[:3],brick,rsdir)

def get_final_dir(outdir,brick,rowstart, 
                  do_skipids='no',do_more='no'):
    """Returns paths like outdir/replaceme/bri/brick/rs0"""
    rsdir= get_rsdir(rowstart, 
                     do_skipids,do_more)
    return os.path.join(outdir,'replaceme',brick[:3],brick,
                        rsdir)


def get_deldirs(outdir,brick,rowstart, 
                do_skipids='no',do_more='no'):
    """If slurm timeout or failed, logfile will exist in final dir but other outputs
        will be in interm dir. Return list of dirst to all of these
    """
    dirs= [get_final_dir(outdir,brick,rowstart, 
                          do_skipids,do_more).replace('replaceme','logs')]
    dirs+= [get_interm_dir(outdir,brick,rowstart, 
                          do_skipids,do_more)]
    return dirs

def get_logdir(outdir,brick,rowstart, 
               do_skipids='no',do_more='no'):
   return (get_final_dir(outdir,brick,rowstart, 
                        do_skipids,do_more)
           .replace('replaceme','logs'))


def get_logfile(outdir,brick,rowstart, 
                do_skipids='no',do_more='no'):
    logdir= get_logdir(outdir,brick,rowstart, 
                       do_skipids,do_more)
    return os.path.join(logdir,'log.'+brick)


def get_slurm_files(outdir):
    return glob( outdir + '/slurm-*.out')
 


class QdoList(object):
    """Queries the qdo db and maps log files to tasks and task status

    Args:
        outdir: obiwan outdir, the slurm*.out files are there
        que_name: ie. qdo create que_name
        skip_suceeded: number succeeded tasks can be very large for production runs, 
            this slows down code so skip those tasks
    """
    def __init__(self,outdir,que_name='obiwan',
                 skip_succeed=False,
                 rand_1e3=False,
                 firstN=None):
        print('que_name= ',que_name.upper())
        self.outdir= outdir
        self.que_name= que_name
        self.skip_succeed= skip_succeed
        self.rand_1e3= rand_1e3
        self.firstN= firstN

    def get_tasks_logs(self):
        """get tasks and logs for the three types of qdo status"""
        # Logs for all Failed tasks
        tasks={}
        ids={}
        logs= defaultdict(list)
        #err= defaultdict(lambda: [])
        q = qdo.connect(self.que_name)
        for res in QDO_RESULT:
            if self.skip_succeed and res == 'succeeded':
                continue
            # List of "brick rs" for each QDO_RESULT 
            qdo_tasks= np.array(q.tasks(state= getattr(qdo.Task, res.upper())))
            print(type(qdo_tasks))
            if self.rand_1e3:
                qdo_tasks= qdo_tasks[np.random.randint(0,len(qdo_tasks),size=1000)]
            elif not self.firstN is None:
                qdo_tasks= qdo_tasks[:self.firstN]
            ids[res],tasks[res] = zip(*[(a.id,a.task) 
                                         for a in qdo_tasks])
            # Corresponding log, slurm files  
            for task in tasks[res]:
                # Logs
                brick,rs,do_skipids,do_more = task.split(' ')
                logfn= get_logfile(self.outdir,brick,rs, 
                                   do_skipids=do_skipids,do_more=do_more)
                logs[res].append( logfn )
        return tasks,ids,logs

    def change_task_state(self,task_ids,to=None, modify=False):
        """change qdo tasks state, for tasks with task_ids, to pending,failed, etc

        Args:
          to: change qdo state to this, pending,failed
          modify: True to actually reset the qdo tasks state AND to delete
            all output files for that task
        """
        assert(to in ['pending','failed'])
        q = qdo.connect(self.que_name)
        for task_id in task_ids:
            try:
                task_obj= q.tasks(id= int(task_id))
                brick,rs,do_skipids,do_more= task_obj.task.split(' ')
                del_dirs= get_deldirs(self.outdir,brick,rs, 
                                       do_skipids=do_skipids,
                                       do_more=do_more)
                if modify:
                    if to == 'pending':
                        # Stuck in pending b/c slurm job timed out
                        task_obj.set_state(qdo.Task.PENDING)
                        # rerun from scratch
                        for dr in del_dirs:
                            dobash('rm -r %s/*' % dr)
                    elif to == 'failed':
                        # Manually force to failed, keep whatever outputs have
                        task_obj.set_state(qdo.Task.FAILED)
                else:
                    print('set --modify to affect id=%d, which corresponds to taks_obj=' % task_id,task_obj)
                    print('if pending, would removed:',del_dirs)
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
            r'ValueError:\ starting\ row=[0-9]*?\ exceeds.*?sources',
            r'\No\ randoms\ in\ brick',
            #'WARING: found nothing with:',
            #'File "obiwan/kenobi.py", line 112, in get_skip_ids',
            #'no skippedids.fits files exist for this brick',
            #"KeyError: 'do_skipds'",
            r'psycopg2\.OperationalError:',
            r'MemoryError',
            r'RuntimeError:\ Command\ failed:\ sex\ -c',
            r'SystemError:\ \<built-in\ method\ flush']
        self.regex_errs_extra= ['Other','log does not exist']

    def get_tally(self):
        tally= defaultdict(list)
        for res in ['succeeded','failed','running']:
            print('res=%s' % res)
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
                for cnt,log in enumerate(self.logs[res]):
                    if (cnt+1) % 25 == 0: print('%d/%d' % (cnt+1,len(self.logs[res])))
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
    parser.add_argument('--skip_succeed',action='store_true',default=False,help='speed up, number succeeded tasks can be very large for production runs and slows down status code',required=False)
    parser.add_argument('--rand_1e3',action='store_true',default=False,help='speed up, look at up to 1000 random succeed,failed,running tasks to get sense of failure modes occuring',required=False)
    parser.add_argument('--firstN',type=int,default=None,help='speed up, instead of random 1000 do the first N (user specified) qdo tasks',required=False)
    parser.add_argument('--running_to_pending',action="store_true",default=False,help='set to reset all "running" jobs to "pending"')
    parser.add_argument('--running_to_failed',action="store_true",default=False,help='set to reset all "running" jobs to "failed"')
    parser.add_argument('--failed_message_to_pending',action='store',default=None,help='set to message of failed tak and reset all failed tasks with that message to pending')
    parser.add_argument('--modify',action='store_true',default=False,help='set to actually reset the qdo tasks state AND to delete IFF running_to_pending or failed_message_to_pending are set')
    parser.add_argument('--outdir',default='.',help='',required=False)
    args = parser.parse_args()
    print(args)

    Q= QdoList(args.outdir,que_name=args.qdo_quename,
               skip_succeed=args.skip_succeed,
               rand_1e3=args.rand_1e3,
               firstN=args.firstN)
    print('Getting tasks,logs')
    tasks,ids,logs= Q.get_tasks_logs()

    # Write log fns so can inspect
    for res in logs.keys():
        writelist(logs[res],"%s_%s_logfns.txt" % (args.qdo_quename,res))
    R= RunStatus(tasks,logs)
    print('Counting number of failed,suceed,running tasks for ech failure mode')
    tally= R.get_tally()
    R.print_tally(tally)

    # logs,tasks for each type of failure
    print('Writing logs,tasks for each failure mode for failed tasks')
    for err_key in R.regex_errs + R.regex_errs_extra:
        err_logs= np.array(logs['failed'])[ tally['failed'] == err_key ]
        err_tasks= np.array(tasks['failed'])[ tally['failed'] == err_key ]
        err_string= (err_key[:12] + err_key[-8:])
        err_string= ((err_key[:10] + err_key[-10:])
                     .replace(" ","_")
                     .replace("/","")
                     .replace("*","")
                     .replace("?","")
                     .replace(":",""))
        writelist(err_logs,"logs_%s_%s.txt" % (args.qdo_quename,err_string))
        writelist(err_tasks,"tasks_%s_%s.txt" % (args.qdo_quename,err_string))

    # Rerun tasks and delete those tasks' outputs
    if len(ids['running']) > 0:
        if args.running_to_pending:
            Q.change_task_state(ids['running'], to='pending',modify=args.modify)
        elif args.running_to_failed:
            Q.change_task_state(ids['running'], to='failed',modify=args.modify)
    
    if args.failed_message_to_pending:
        hasMessage= np.where(tally['failed'] == args.failed_message_to_pending)[0]
        if hasMessage.size > 0:
            theIds= np.array(ids['failed'])[hasMessage]
            Q.change_task_state(theIds, to='pending', modify=args.modify)

    print('done')
