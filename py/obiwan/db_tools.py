'''
Query our PostgreSQL Databases at NERSC
'''

import psycopg2
import os
import numpy as np

try:
    from astrometry.util.fits import fits_table, merge_tables
except ImportError:
    pass

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

def all_psqlcols_for_ids(ids, db_randoms_table='obiwan_elg_ra175',
                       try_with_join=False):
    """Returns all db columns in the db having the ids provided

    Args:
        ids: list or array, ids generally come from obiwan 'simcat*.fits' table, for example
        db_table: table name in psql db 'desi' hosted at 'scidb2.nersc.gov'
        try_with_join: to use equivalent sql select that uses join 
    """
    db= PsqlWorker()
    columns= 'id ra dec g r z rhalf n ba pa redshift'.split(' ')
    if db_randoms_table == 'eboss_elg':
        columns+= 'id_sample,nn_redshift'.split(' ')
    cmd= "SELECT "
    for col in columns[:-1]:
        cmd+= "%s," % col
    cmd+= "%s" % columns[-1]

    if not try_with_join: 
        # Simplest, faster
        cmd+= " FROM %s WHERE id in (" % db_randoms_table
        for i in ids[:-1]:
            cmd+= "%d," % i
        cmd+= "%d)" % ids[-1]
    else:
        # Slower
        vals=""
        for i in ids[:-1]:
            vals+= "(%d)," % i
        vals+= "(%d)" % ids[-1]
        cmd= (cmd + " FROM %s as db RIGHT JOIN (values %s) " 
              % (db_randoms_table,vals) 
              + 
              "as v(id) on (db.id=v.id)")
    print('cmd= %s' % cmd)
    db.cur.execute(cmd)
    # List of tuples [(id,reshift,...),(id,reshift,...)]
    a= db.cur.fetchall() 
    # Tuple of lists (ids,reshifts,...)
    tup= zip(*a)
    #tup[ith_col])
    return {col: np.array(vals) 
            for col,vals in zip(columns,tup)}
    #return np.array(sql_ids),np.array(sql_redshift)

if __name__ == '__main__':
    T= getSrcsInBrick('1765p247','elg', db_table='obiwan_elg_ra175')

    simcat= fits_table("/global/cscratch1/sd/kaylanb/obiwan_out/elg_9deg2_ra175/elg/176/1765p247/rs0/obiwan/simcat-elg-1765p247.fits")
    data_dict= all_psqlcols_for_ids(simcat.id, db_randoms_table='obiwan_elg_ra175')
    for i in range(10):
        print(data_dict['id'][i],data_dict['redshift'][i])
    
