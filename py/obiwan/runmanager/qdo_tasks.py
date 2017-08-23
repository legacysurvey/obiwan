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
    I,J,d = match_radec(brick.ra, brick.dec, box_center['ra'],box_center['dec'],
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
  b= Bricks()
  tab= b.overlapBox(ra=[123.3,124.3], dec=[24.0,25.0])

  from obiwan.runmanager.status import writelist
  writelist(tab.brickname, 'bricks_10deg.txt')
