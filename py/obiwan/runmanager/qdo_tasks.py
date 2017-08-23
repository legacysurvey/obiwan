import os
import numpy as np

from astrometry.util.fits import fits_table
from astrometry.libkd.spherematch import match_radec

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

if __name__ == '__main__':
  from argparse import ArgumentParser
  parser = ArgumentParser()
  parser.add_argument('--ra1',type=float,action='store',default=123.3,required=False)
  parser.add_argument('--ra2',type=float,action='store',default=124.3,required=False)
  parser.add_argument('--dec1',type=float,action='store',default=24.0,required=False)
  parser.add_argument('--dec2',type=float,action='store',default=25.0,required=False)
  parser.add_argument('--nobj_total',type=int,action='store',default=2400,help='total number of randoms in the region running obiwan on',required=False)
  parser.add_argument('--nobj_per_run',type=int,action='store',default=500,help='number of simulated sources to inject per obiwan run',required=False)
  args = parser.parse_args()

  # bricks touching region
  b= Bricks()
  tab= b.overlapBox(ra=[args.ra1,args.ra2], dec=[args.dec1,args.dec2])

  from obiwan.runmanager.status import writelist
  writelist(np.sort(tab.brickname), 'bricks_inregion.txt')

  # corresponding qdo tasks
  avg_nobj_per_brick= int( args.nobj_total / len(tab))
  print('avg_nobj_per_brick= %d' % avg_nobj_per_brick)
  tasks= ['%s %d' % (brick,rs) 
            for brick in np.sort(tab.brickname)
            for rs in np.arange(0,avg_nobj_per_brick,args.nobj_per_run)]
  writelist(tasks, 'tasks_inregion.txt')
