import os
import numpy as np

from astrometry.util.fits import fits_table
from astrometry.libkd.spherematch import match_radec

from obiwan.runmanager.status import writelist

SURVEY_BRICKS= os.path.join(os.environ['obiwan_data'],
                            'legacysurveydir/survey-bricks.fits.gz')
assert(os.path.exists(SURVEY_BRICKS))

class Bricks(object):
    def __init__(self):
        self.bricks= fits_table(SURVEY_BRICKS)

    def overlapBox(self,ra=[100,101],dec=[20,21]):
        """ra,dec: corners of box injected randoms into"""
        hw=0.25
        return self.bricks[((ra[0] -hw <= self.bricks.ra) &
                            (ra[1] +hw >= self.bricks.ra) &
                            (dec[0] -hw <= self.bricks.dec) &
                            (dec[1] +hw >= self.bricks.dec))

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

    def get_bricks(self):
        """Returns bricks in ra,dec region"""
        b= Bricks()
        self.bricks= b.overlapBox(ra=[self.ra1,self.ra2], dec=[self.dec1,self.dec2])
        #writelist(np.sort(self.bricks.brickname), 'bricks_inregion.txt')

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

    def tasklist_skipids(self,do_more='no',minid=None):
        """tasklist for skipids runs

        Args:
            do_more: yes or no, yes if running more randoms b/c TS returns too few target
            minid: if do_more == yes, this must be an integer for the randoms id to start from
        """
        do_skipids= 'yes'
        tasks= [self.task(brick,rs,do_skipid,do_more) 
                for brick in np.sort(self.bricks.brickname)
                for rs in np.arange(0,3*self.nobj_per_run, self.nobj_per_run)]
        writelist(tasks, 'tasks_skipid_%s_more_%s_minid_%s.txt' % 
                         (do_skipids,do_more,str(minid)))


    def get_tasklist(self,objtype='elg',randoms_db='obiwan_elg_ra175',
                     do_more='no',minid=1,outdir=None,
                     bricks=None,estim_nperbrick=2e3):
        """It is too slow to find the example number of randoms in each brick, so find all bricks
            with at least one source and put in the expected number of randoms + 2 StdErros worth

        Args:
            objtype: elg,lrg
            minid: None, unless do_more == yes then it is an integer for the randoms id to start from
            bricks: array like list of bricks to get qdo tasks for, if None all bricks found
        """
        tasks=[]
        if bricks is None:
            bricks= np.sort(self.bricks.brickname)
        else:
            bricks= np.sort(bricks)
        tasks+= [self.task(brick,rs,do_skipids,do_more) 
                 for rs in np.arange(0,estim_nperbrick,self.nobj_per_run)]
        return tasks

    def writetasks(self,tasks,not_in_bricks,
                   do_more='no',minid=1,do_skipids='no'):
        """Write task list to file"""
        writelist(tasks, 'tasks_skipid_%s_more_%s_minid_%s.txt' % 
                    (do_skipids,do_more,str(minid)))
        print('Number of bricks without randoms: %d' % len(not_in_bricks))
        print('First 10 are',not_in_bricks[:10])


if __name__ == '__main__':
    do_skipids='no'
    do_more='yes'
    if do_more == 'yes':
        minid=240001
    else:
        minid=None
    ###
    d= dict(ra1=150.,ra2=160.,
            dec1=0.,dec2=10.0,
            nobj_total=2.4e6,nobj_per_run=300)
    outdir='/global/cscratch1/sd/kaylanb/obiwan_out/elg_100deg2'
    objtype='elg'
    # Initialize
    T= TaskList(**d)
    T.get_bricks()
    num= T.estim_nperbrick()
    # Write tasks
    if do_skipids == 'no':
        T.tasklist(objtype,
                   do_more,minid,outdir,
                   estim_nperbrick=num)
    else:
        T.tasklist_skipids(do_more=do_more,minid=minid)
    # Qdo "skip_ids" task list for a given list of bricks
    #brick_list_fn= '/global/cscratch1/sd/kaylanb/obiwan_code/obiwan/bricks_ready_skip.txt'
    #write_qdo_tasks_skipids(brick_list_fn, nobj_per_run=300)
    #from argparse import ArgumentParser
    #parser = ArgumentParser()
    #parser.add_argument('--region', type=str, choices=['no','yes'],default='no', help='inject skipped ids for brick, otherwise run as usual')
    #parser.add_argument('--bricks_list', type=str, default=None, help='text file listing bricks to inject skipped ids for brick, otherwise run as usual')
    #parser.add_argument('--ra1',type=float,action='store',default=123.3,required=False)
    #parser.add_argument('--ra2',type=float,action='store',default=124.3,required=False)
    #parser.add_argument('--dec1',type=float,action='store',default=24.0,required=False)
    #parser.add_argument('--dec2',type=float,action='store',default=25.0,required=False)
    #parser.add_argument('--nobj_total',type=int,action='store',default=2400,help='total number of randoms in the region running obiwan on',required=False)
    #parser.add_argument('--nobj_per_run',type=int,action='store',default=500,help='number of simulated sources to inject per obiwan run',required=False)
#args = parser.parse_args()

