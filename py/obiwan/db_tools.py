'''
function to query PSQL db for all sources in a Brick
'''

import psycopg2
import os
import numpy as np

from astrometry.util.fits import fits_table, merge_tables

class PsqlWorker(object):
    def __init__(self):
        self.conn= psycopg2.connect(host='nerscdb03.nersc.gov', user='desi_admin', database='desi')
        self.cur= self.conn.cursor()

    def close(self):
        self.cur.close()
        self.conn.close()

def getSrcsInBrick(brickname,objtype, db_table='obiwan_elg',
                   skipped_ids=None):
    """Returns tuple: fits table, seed
    
    Args:
        skipped_ids: array or list of strings of ids if not None, the db ids
    """
    db= PsqlWorker()
    cmd= "select brickname,brickid,ra1,ra2,dec1,dec2 from obiwan_bricks where brickname = '%s'" % brickname
    print('cmd= %s' % cmd)
    db.cur.execute(cmd)
    a= db.cur.fetchall()
    assert(len(a) == 1)
    name,brickid,ra1,ra2,dec1,dec2= a[0]
    assert(name == brickname)
    #ra1,ra2,dec1,dec2= 10,50,10,20
    cmd="select id,ra,dec,g,r,z"
    if objtype in ['elg','lrg']:
        cmd=cmd+ ",redshift,rhalf,n,ba,pa"
        if objtype == 'lrg':
            cmd=cmd+ ",w1" % objtype
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
        raise ValueError('No randoms in brick %s, e.g. found nothing with db query: %s' % (brickname,cmd))
    # Package in fits_table
    d={}
    # TODO: make simpler and use re instead of rhalf above
    if objtype == 'star':
        d['id'],d['ra'],d['dec'],d['g'],d['r'],d['z']= zip(*a)
    elif objtype == 'elg':
        d['id'],d['ra'],d['dec'],d['g'],d['r'],d['z'],\
                d['redshift'],d['rhalf'],d['n'],d['ba'],\
                d['pa']= zip(*a)
    elif objtype == 'lrg':
        d['id'],d['ra'],d['dec'],d['g'],d['r'],d['z'],\
                d['redshift'],d['rhalf'],d['n'],d['ba'],\
                d['pa'],d['w1']= zip(*a)
    del a
    T= fits_table()
    for key in d.keys():
        T.set(key, np.array(d[key]))
    del d
    #
    db.close()
    return T,brickid

def redshifts_for_ids(ids, db_table='obiwan_elg_ra175',
                     try_with_join=False):
    """Returns the reshifts of randoms in the db having the ids provided

    Args:
        ids: list or array, ids generally come from obiwan 'simcat*.fits' table, for example
        db_table: table name in psql db 'desi' hosted at 'scidb2.nersc.gov'
        try_with_join: to use equivalent sql select that uses join 
    """
    db= PsqlWorker()
    # Simplest
    cmd= "SELECT id,elg_redshift FROM %s WHERE id in (" % db_table
    for i in ids[:-1]:
        cmd+= "%d," % i
    cmd+= "%d)" % ids[-1]
    if try_with_join:
        # Equvialent but diff speed 
        vals=""
        for i in ids[:-1]:
            vals+= "(%d)," % i
        vals+= "(%d)" % ids[-1]
        cmd= ("SELECT db.id,db.elg_redshift FROM %s as db RIGHT JOIN (values %s) " 
              % (db_table,vals) 
              + 
              "as v(id) on (db.id=v.id)")
    print('cmd= %s' % cmd)
    db.cur.execute(cmd)
    a= db.cur.fetchall() #list of tuples (id,reshift)
    #a= zip(*a)
    #sql_ids,sql_redshift= a[0],a[1]
    sql_ids,sql_redshift= zip(*a)
    return np.array(sql_ids),np.array(sql_redshift)

if __name__ == '__main__':
    T= getSrcsInBrick('1765p247','elg', db_table='obiwan_elg_ra175')

    simcat= fits_table("/global/cscratch1/sd/kaylanb/obiwan_out/elg_9deg2_ra175/elg/176/1765p247/rs0/obiwan/simcat-elg-1765p247.fits")
    ids,redshifts= redshifts_for_ids(simcat.id, db_table='obiwan_elg_ra175')
    for id,z in zip(ids[:10],redshifts[:10]):
        print(id,z)
    
