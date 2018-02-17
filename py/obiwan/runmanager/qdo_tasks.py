import os
import numpy as np

from astrometry.util.fits import fits_table

from obiwan.common import writelist


class Bricks(object):
    def __init__(self,survey_bricks_fn):
        self.bricks= fits_table(survey_bricks_fn)

    def overlapBox(self,ra=[100,101],dec=[20,21]):
        """ra,dec: corners of box injected randoms into"""
        hw=0.25
        return self.bricks[((ra[0] -hw <= self.bricks.ra) &
                            (ra[1] +hw >= self.bricks.ra) &
                            (dec[0] -hw <= self.bricks.dec) &
                            (dec[1] +hw >= self.bricks.dec))]

class TaskList(object):
    """Creates QDO tasks lists for default,do_skipids, and/or do_more

    Args:
        ra1,ra2,dec1,dec2: floats, corners of ra dec box
        nobj_per_run: number of simulated sources to inject per obiwan run
    """

    def __init__(self,ra1=123.3,ra2=124.3,dec1=24.0,dec2=25.0,
                 nobj_total=1e6,nobj_per_run=500): 
        self.ra1=ra1
        self.ra2=ra2
        self.dec1=dec1
        self.dec2=dec2
        self.nobj_total=nobj_total
        self.nobj_per_run=nobj_per_run

    def bricks_in_region(self,survey_bricks_fn=None):
        """Returns bricks in ra,dec region"""
        b= Bricks(survey_bricks_fn=survey_bricks_fn)
        self.bricks= b.overlapBox(ra=[self.ra1,self.ra2], dec=[self.dec1,self.dec2])
        #writelist(np.sort(self.bricks.brickname), 'bricks_inregion.txt')

    def bricks_from_file(self,fn):
        self.bricks= np.loadtxt(fn,dtype=str)

    def estim_nperbrick(self):
        p= 1./len(self.bricks)
        Ex= self.nobj_total * p
        SE= np.sqrt(self.nobj_total * p * (1-p))
        return Ex + 2*SE

    def task(self,brick,rs,do_skipids,do_more):
        """returns a single QDO task as a string"""
        assert(do_skipids in ['yes','no'])
        assert(do_more in ['yes','no'])
        return '%s %d %s %s' % (brick,rs,do_skipids,do_more) 

    def tasklist_skipids(self,bricks=None,do_more='no',minid=None):
        """tasklist for skipids runs

        Args:
            do_more: yes or no, yes if running more randoms b/c TS returns too few target
            minid: if do_more == yes, this must be an integer for the randoms id to start from
        """
        do_skipids= 'yes'
        if bricks is None:
            bricks= np.sort(self.bricks.brickname)
        tasks= [self.task(brick,rs,do_skipid,do_more) 
                for brick in bricks
                for rs in np.arange(0,3*self.nobj_per_run, self.nobj_per_run)]
        writelist(tasks, 'tasks_skipid_%s_more_%s_minid_%s.txt' % 
                         (do_skipids,do_more,str(minid)))


    def get_tasklist(self,bricks=None,objtype='elg',
                     do_more='no',minid=1,
                     estim_nperbrick=2e3,
                     cosmos=False):
        """It is too slow to find the example number of randoms in each brick, so find all bricks
            with at least one source and put in the expected number of randoms + 2 StdErros worth

        Args:
            objtype: elg,lrg
            minid: None, unless do_more == yes then it is an integer for the randoms id to start from
            bricks: array like list of bricks to get qdo tasks for, if None all bricks found
        """
        #tasks=[]
        if bricks is None:
            bricks= np.sort(self.bricks.brickname)
        else:
            bricks= np.sort(bricks)
        tasks= [self.task(brick,rs,do_skipids,do_more) 
                for rs in np.arange(0,estim_nperbrick,self.nobj_per_run)
                for brick in bricks]
        return tasks

    def writetasks(self,tasks,
                   do_more='no',minid=1,do_skipids='no'):
        """Write task list to file"""
        fn= 'tasks_skipid_%s_more_%s_minid_%s.txt' % \
                (do_skipids,do_more,str(minid))
        writelist(tasks, fn)


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--obj', type=str,choices=['elg','star'],required=True)
    parser.add_argument('--radec', nargs='+',type=float, 
                        help='no quotes, e.g. --radec ra1 ra2 dec1 dec2',required=True)
    parser.add_argument('--nobj_total', type=int, default=10000,required=True)
    parser.add_argument('--survey_bricks_fn', type=str, 
                        help='abs path to survey-bricks.fits.fz',
                        required=True)
    parser.add_argument('--nobj_per_run', type=int, default=300)
    parser.add_argument('--bricks_fn', type=str, default=None,
                        help='abs path to list of bricks for the tasklist')
    parser.add_argument('--do_skipids', type=str, choices=['yes','no'],
                        default='no')
    parser.add_argument('--do_more', type=str, choices=['yes','no'],
                        default='no')
    parser.add_argument('--minid', type=int, default=None)
    parser.add_argument('--use_bricklist_given', action='store_true', default=False)
    parser.add_argument('--cosmos', action='store_true', default=False,
                        help='set to add the cosmos subset number to the tasks')
    args = parser.parse_args()

    assert(len(args.radec) == 4)
    if args.do_more == 'yes':
        assert(not args.minid is None)
    do_skipids='no'
    do_more='no'
    #survey_bricks=os.path.join(os.environ['HOME'],
    #                        'mydata/survey-bricks.fits.gz')
    #survey_bricks=os.path.join(os.environ['HOME'],
    #                        'Downloads/survey-bricks-dr5.fits.gz')
    #survey_bricks= os.path.join(os.environ['obiwan_data'],
    #                        'legacysurveydir/survey-bricks.fits.gz')
    
    ###
    d=dict(nobj_total=args.nobj_total,
           nobj_per_run=args.nobj_per_run)
    if not args.use_bricklist_given:
        d.update(ra1=args.radec[0],ra2=args.radec[1],
                 dec1=args.radec[2],dec2=args.radec[3])
    # Initialize
    T= TaskList(**d)
    if args.use_bricklist_given:
        T.bricks_from_file(args.bricks_fn)
    else:
        T.bricks_in_region(survey_bricks_fn=args.survey_bricks_fn)
    num= T.estim_nperbrick()
    # Write tasks
    if args.do_skipids == 'no':
        tasks= T.get_tasklist(T.bricks,args.obj,
                              args.do_more,args.minid,
                              estim_nperbrick=num)
    else:
        T.tasklist_skipids(T.bricks,
                           do_more=do_more,minid=minid)
    if args.cosmos:
        tasks= ["%s %d" % (task,subset)
                for task in tasks
                for subset in [60,64,69]]
    T.writetasks(tasks,
                 do_more=args.do_more,minid=args.minid,
                 do_skipids=args.do_skipids)
    # Qdo "skip_ids" task list for a given list of bricks
    #brick_list_fn= '/global/cscratch1/sd/kaylanb/obiwan_code/obiwan/bricks_ready_skip.txt'
    #write_qdo_tasks_skipids(brick_list_fn, nobj_per_run=300)
    

