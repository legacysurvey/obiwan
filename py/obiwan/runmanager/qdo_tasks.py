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
    """approx: bricks within radius of hypotenjuse of box of box center"""
    box_center=dict(ra= np.average(ra),
                    dec= np.average(dec))
    rad_deg= np.sqrt((box_center['ra'] - ra[0])**2 + \
                     (box_center['dec'] - dec[0])**2)
    print('box_center=',box_center,'rad_deg=',rad_deg)
    I,J,d = match_radec(self.bricks.ra, self.bricks.dec, 
                        box_center['ra'],box_center['dec'],
                        rad_deg, nearest=True)
    return self.bricks[I].copy()
    #for corn_ra,corn_dec in zip([(ra[0],dec[0]),
    #                             (ra[1],dec[0]),
    #                             (ra[0],dec[1]),
    #                             (ra[1],dec[1])]):
    #keep= np.logical_or.reduce((self.bricks.ra1 >= ra[0],
    #                             self.bricks.ra2 <= ra[1],
    #                             self.bricks.dec1 >= dec[0],
    #                             self.bricks.dec2 <= dec[1]))
    #return self.bricks[keep].copy()

class TaskList(object):
    """Creates QDO tasks lists for default,do_skipids, and/or do_more

    Args:
        ra1,ra2,dec1,dec2: floats, corners of ra dec box
        nobj_per_run: number of simulated sources to inject per obiwan run
    """

    def __init__(self,ra1=123.3,ra2=124.3,dec1=24.0,dec2=25.0,
                 nobj_per_run=500): 
        self.ra1=ra1
        self.ra2=ra2
        self.dec1=dec1
        self.dec2=dec2
        self.nobj_per_run=nobj_per_run

    def get_bricks(self):
        """Returns bricks in ra,dec region"""
        b= Bricks()
        self.bricks= b.overlapBox(ra=[self.ra1,self.ra2], dec=[self.dec1,self.dec2])
        writelist(np.sort(self.bricks.brickname), 'bricks_inregion.txt')

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

    def tasklist(self,objtype,randoms_db,
                 do_more,minid,outdir):
        """For each brick, gets all randoms from PSQL db, and finds exact number of rs* tasks needed

        Args:
            objtype: elg,lrg
            randoms_db: name of PSQL db for randoms, e.g. obiwan_elg_ra175
            do_more: yes or no, yes if running more randoms b/c TS returns too few target
            minid: None, unless do_more == yes then it is an integer for the randoms id to start from
            outdir: path/to/obiwan_out/
        """
        do_skipids='no'
        from obiwan.kenobi import get_sample
        sample_kwargs= {"objtype":objtype,
                        "randoms_db":randoms_db,
                        "minid":minid,
                        "outdir":outdir,
                        "do_skipids":do_skipids}
        tasks=[]
        bricks= np.sort(self.bricks.brickname)
        not_in_bricks=[]
        for cnt,brick in enumerate(bricks):
            if cnt % 10 == 0: print('brick %d/%d' % (cnt+1,len(bricks)))
            try: 
                Samp= get_sample(brick=brick,verbose=False,**sample_kwargs)
            except ValueError as err_obj:
                if 'No randoms in brick' in str(err_obj):
                    not_in_bricks+= [brick]
                    continue
                else: 
                    raise ValueError(str(err_obj))
            tasks+= [self.task(brick,rs,do_skipids,do_more) 
                     for rs in np.arange(0,len(Samp),self.nobj_per_run)]
        writelist(tasks, 'tasks_skipid_%s_more_%s_minid_%s.txt' % 
                         (do_skipids,do_more,str(minid)))
        print('number not in bricks = %d' % len(not_in_bricks))
        print(not_in_bricks[:30])


if __name__ == '__main__':
    do_skipids='no'
    do_more='yes'
    if do_more == 'yes':
    minid=240001
    else:
    minid=None
    ###
    ra1=173.5
    ra2=176.5
    dec1=23.0
    dec2=26.0
    randoms_db='obiwan_elg_ra175'
    outdir='/global/cscratch1/sd/kaylanb/obiwan_out/elg_9deg2_ra175'
    nobj_per_run=300
    objtype='elg'
    # Initialize
    T= TaskList(ra1=ra1,ra2=ra2, dec1=dec1,dec2=dec2,
              nobj_per_run=nobj_per_run)
    T.get_bricks()
    # Write tasks
    if do_skipids == 'no':
        T.tasklist(objtype,randoms_db,
                   do_more,minid,outdir)
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

