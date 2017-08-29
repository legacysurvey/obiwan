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

def write_qdo_tasks_normal(ra1=123.3,ra2=124.3,dec1=24.0,dec2=25.0, 
                           nobj_total=2400, nobj_per_run=500):
  """for all bricks in a ra,dec box region, write qdo task list for obiwan runs
  
  Args:
    ra1,ra2,dec1,dec2: floats, corners of ra dec box
    nobj_total: total number of randoms in the region running obiwan on
    nobj_per_run: number of simulated sources to inject per obiwan run

  Returns:
    Writes qdo task list to file
  """
  # bricks touching region
  b= Bricks()
  tab= b.overlapBox(ra=[ra1,ra2], dec=[dec1,dec2])

  writelist(np.sort(tab.brickname), 'bricks_inregion.txt')

  # corresponding qdo tasks
  avg_nobj_per_brick= int( nobj_total / len(tab))
  print('avg_nobj_per_brick= %d' % avg_nobj_per_brick)
  # Tasks: brickname rowstart skip_ids
  tasks= ['%s %d %s' % (brick,rs,'no') 
            for brick in np.sort(tab.brickname)
            for rs in np.arange(0,avg_nobj_per_brick, nobj_per_run)]
  writelist(tasks, 'tasks_inregion.txt')

def write_qdo_tasks_skipids(brick_list_fn, nobj_per_run=500):
  """for given list of bricks, write qdo task list for "skip_ids" obiwan runs

  Args:
    brick_list_fn: text file with one brickname per line
    nobj_per_run: number of simulated sources to inject per obiwan run

  Returns:
    Writes qdo task list to file
  """
  bricks= np.loadtxt(brick_list_fn,dtype=str)
  tasks= ['%s %d %s' % (brick,rs,'yes') 
            for brick in bricks
            for rs in np.arange(0,2*nobj_per_run, nobj_per_run)]
  writelist(tasks, 'tasks_skipids.txt')



if __name__ == '__main__':
  # Qdo "normal" task list for all bricks in a ra,dec box region
  ##write_qdo_tasks_normal(ra1=123.3,ra2=124.3,dec1=24.0,dec2=25.0, 
  ##                       nobj_total=2400, nobj_per_run=500)
  # Qdo "skip_ids" task list for a given list of bricks
  brick_list_fn= '/global/cscratch1/sd/kaylanb/obiwan_code/obiwan/bricks_ready_skip.txt'
  write_qdo_tasks_skipids(brick_list_fn, nobj_per_run=300)
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

