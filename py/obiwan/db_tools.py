'''
function to query PSQL db for all sources in a Brick
'''

import psycopg2
import os
import numpy as np
from astrometry.util.fits import fits_table, merge_tables

class PsqlWorker(object):
    def __init__(self):
        self.conn= psycopg2.connect(host='scidb2.nersc.gov', user='desi_admin', database='desi')
        self.cur= self.conn.cursor()

    def close(self):
        self.cur.close()
        self.conn.close()

def getSrcsInBrick(brickname,objtype, db_table='obiwan_elg',
                   skipped_ids=None):
    '''returns a fit_table of PSQL query
    
    Args:
      skipped_ids: array or list of strings of ids if not None, the db ids
    '''
    db= PsqlWorker()
    cmd= "select brickname,ra1,ra2,dec1,dec2 from obiwan_bricks where brickname = '%s'" % brickname
    print('cmd= %s' % cmd)
    db.cur.execute(cmd)
    a= db.cur.fetchall()
    assert(len(a) == 1)
    name,ra1,ra2,dec1,dec2= a[0]
    assert(name == brickname)
    #ra1,ra2,dec1,dec2= 10,50,10,20
    cmd="select id,seed,ra,dec,%s_g,%s_r,%s_z" % \
        tuple([objtype]*3)
    if objtype in ['elg','lrg']:
        cmd=cmd+ ",%s_redshift,%s_rhalf,%s_n,%s_ba,%s_pa" % \
                 tuple([objtype]*5)
        if objtype == 'lrg':
            cmd=cmd+ ",%s_w1" % objtype
    cmd=cmd+ " from %s where q3c_poly_query(ra, dec, '{%f,%f, %f,%f, %f,%f, %f,%f}')" % \
             (db_table, ra1,dec1, ra2,dec1, ra2,dec2, ra1,dec2)
    if skipped_ids is not None:
      cmd+= "and id in ("
      for skip_id in skipped_ids[:-1]:
        cmd+= "%s," % skip_id
      cmd+= "%s)" % skipped_ids[-1]
    print('cmd= %s' % cmd)
    db.cur.execute(cmd)
    a= db.cur.fetchall()
    if len(a) == 0:
        print('WARING: found nothing with: %s' % cmd)
    # Package in fits_table
    d={}
    # TODO: make simpler and use re instead of rhalf above
    if objtype == 'star':
        d['id'],d['seed'],d['ra'],d['dec'],d['star_g'],d['star_r'],d['star_z']= zip(*a)
    elif objtype == 'elg':
        d['id'],d['seed'],d['ra'],d['dec'],d['elg_g'],d['elg_r'],d['elg_z'],\
                d['elg_redshift'],d['elg_re'],d['elg_n'],d['elg_ba'],\
                d['elg_pa']= zip(*a)
    elif objtype == 'lrg':
        d['id'],d['seed'],d['ra'],d['dec'],d['lrg_g'],d['lrg_r'],d['lrg_z'],\
                d['lrg_redshift'],d['lrg_re'],d['lrg_n'],d['lrg_ba'],\
                d['lrg_pa'],d['lrg_w1']= zip(*a)
    del a
    T= fits_table()
    for key in d.keys():
        T.set(key, np.array(d[key]))
    del d
    #
    db.close()
    return T

if __name__ == '__main__':
    #T= getSrcsInBrick('0100p100')
    T= getSrcsInBrick('1238p245')
    
