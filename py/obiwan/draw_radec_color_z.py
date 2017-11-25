# See LICENSE.rst for BSD 3-clause license info
# -*- coding: utf-8 -*-
"""
=========================
obiwan.draw_radec_color_z
=========================

mpi4py to draw N random ra,dec with grzW1,Re,redshift info from KDEs
These N ra,dec rows are written to N/n_tasks fits files
fits2db to load those N/n_tasks fits files into the PostgresQL DB
Add bricks table to DB and index on that
func(brick) -- returns all ra,dec in a given brick
Write n_bricks fits files containing these various ra,dec
"""

from __future__ import division, print_function

import os
import argparse
import numpy as np
from astrometry.util.ttime import Time
import datetime
import sys
from scipy import spatial

from astrometry.util.fits import fits_table, merge_tables

from obiwan.common import fits2pandas

class GaussianMixtureModel(object):
    """
    John's class to read, write, and sample from a mixture model.
    
    Args: 
        weights,means,covars: array-like, from mixture fit
        covar_type: usually 'full'
        py: one of ['27','36']
    """
    def __init__(self, weights, means, covars, 
                 py=None,covar_type='full',is1D=False):
        assert(py in ['27','36'])
        self.py= py
        self.is1D= is1D
        self.weights_ = weights
        self.means_ = means
        self.covariances_ = covars
        if self.is1D:
            self.weights_ = self.weights_.reshape(-1,1)
            self.means_ = self.means_.reshape(-1,1)
            self.covariances_ = self.covariances_.reshape(-1,1,1)
        #    self.n_components, self.n_dimensions = self.means_.shape[0],1
        #else:
        self.n_components, self.n_dimensions = self.means_.shape
        #print(self.weights_.shape,self.covariances_.shape,len(self.covariances_.shape))
        self.covariance_type= covar_type
    
    @staticmethod
    def save(model, filename):
        for name,data in zip(['means','weights','covars'],
                     [model.means_, model.weights_,
                      model.covariances_]):
            fn= '%s_%s.txt' % (filename,name)
            np.savetxt(fn,data,delimiter=',')
            print('Wrote %s' % fn)

    @staticmethod
    def load(name,py=None,is1D=False,indir='./'):
        """name: prefix to _weights.txt or _means.txt"""
        d={key:np.loadtxt(os.path.join(indir,name+'_%s.txt' % key),
                           delimiter=',')
           for key in ['means','weights','covars']}
        return GaussianMixtureModel(
                    d['weights'],d['means'],d['covars'],
                    covar_type='full',py=py,is1D=is1D)
    
    def sample(self, n_samples=1, random_state=None):
        assert(n_samples >= 1)
        self.n_samples= n_samples
        if random_state is None:
            random_state = np.random.RandomState()
        self.rng= random_state
        
        if self.py == '2.7':
            X= self.sample_py2()
        else:
            X,Y= self.sample_py3()
        return X
    
    def sample_py2(self):
        weight_cdf = np.cumsum(self.weights_)
        X = np.empty((self.n_samples, self.n_components))
        rand = self.rng.rand(self.n_samples)
        # decide which component to use for each sample
        comps = weight_cdf.searchsorted(rand)
        # for each component, generate all needed samples
        for comp in range(self.n_components):
            # occurrences of current component in X
            comp_in_X = (comp == comps)
            # number of those occurrences
            num_comp_in_X = comp_in_X.sum()
            if num_comp_in_X > 0:
                X[comp_in_X] = self.rng.multivariate_normal(
                    self.means_[comp], self.covariances_[comp], num_comp_in_X)
        return X
    
    def sample_py3(self):
        """Copied from sklearn's mixture.GaussianMixture().sample()"""
        print(self.weights_.shape)
        try:
            n_samples_comp = self.rng.multinomial(self.n_samples, self.weights_)
        except ValueError:
            self.weights_= np.reshape(self.weights_,len(self.weights_))
            n_samples_comp = self.rng.multinomial(self.n_samples, self.weights_)
        if self.covariance_type == 'full':
            X = np.vstack([
                self.rng.multivariate_normal(mean, covariance, int(sample))
                for (mean, covariance, sample) in zip(
                    self.means_, self.covariances_, n_samples_comp)])
        elif self.covariance_type == "tied":
            X = np.vstack([
                self.rng.multivariate_normal(mean, self.covariances_, int(sample))
                for (mean, sample) in zip(
                    self.means_, n_samples_comp)])
        else:
            X = np.vstack([
                mean + self.rng.randn(sample, n_features) * np.sqrt(covariance)
                for (mean, covariance, sample) in zip(
                    self.means_, self.covariances_, n_samples_comp)])

        y = np.concatenate([j * np.ones(sample, dtype=int)
                           for j, sample in enumerate(n_samples_comp)])

        return (X, y)

def ptime(text,t0):
    tnow=Time()
    print('TIMING:%s ' % text,tnow-t0)
    return tnow

def get_area(radec):
    '''returns area on sphere between ra1,ra2,dec2,dec1
    https://github.com/desihub/imaginglss/model/brick.py#L64, self.area=...
    '''
    deg = np.pi / 180.
    # Wrap around
    if radec['ra2'] < radec['ra1']:
        ra2=radec['ra2']+360.
    else:
        ra2=radec['ra2']
    
    area= (np.sin(radec['dec2']*deg)- np.sin(radec['dec1']*deg)) * \
          (ra2 - radec['ra1']) * \
          deg* 129600 / np.pi / (4*np.pi)
    approx_area= (radec['dec2']-radec['dec1'])*(ra2-radec['ra1'])
    print('approx area=%.2f deg2, actual area=%.2f deg2' % (approx_area,area))
    return area

def get_radec(radec,\
              ndraws=1,random_state=np.random.RandomState()):
    """Draws ndraws samples of Ra,Dec from the unit sphere.

	Args:
		radec: dict with keys ra1,ra2,dec1,dec2
			the ra,dec limits for the sample
		ndraws: number of samples
		randome_state: numpy random number generator

	Returns:
		ra,dec: tuple of arrays having length ndraws

	Note: 
		Taken from 
		https://github.com/desihub/imaginglss/blob/master/scripts/imglss-mpi-make-random.py#L55
	"""
    ramin,ramax= radec['ra1'],radec['ra2']
    dcmin,dcmax= radec['dec1'],radec['dec2']
    u1,u2= random_state.uniform(size=(2, ndraws) )
    #
    cmin = np.sin(dcmin*np.pi/180)
    cmax = np.sin(dcmax*np.pi/180)
    #
    RA   = ramin + u1*(ramax-ramin)
    DEC  = 90-np.arccos(cmin+u2*(cmax-cmin))*180./np.pi
    return RA,DEC

def get_sample_dir(outdir,obj):
    return outdir

def mkdir_needed(d):
    """make each needed directory 
    d= dictionary, vars(args)
    """
    dirs=[d['outdir']]
    dirs.append( get_sample_dir(d['outdir'],d['obj']) )
    for dr in dirs:
        if not os.path.exists(dr):
            os.makedirs(dr)

def write_calling_seq(d):
    """each `*_randoms/` directory should have a file listing how randoms were created

    Args:
      d: dict, vars(args)
    """
    dr= get_sample_dir(d['outdir'],d['obj'])
    fn=os.path.join(dr,'README.txt')
    if os.path.exists(fn):
        os.remove(fn)
    with open(fn,'w') as foo:
        for key in d.keys():
            foo.write('%s %s\n' % (key,str(d[key])))
    print('Wrote %s' % fn)


def get_sample_fn(seed=None):
    return 'randoms_rank_%d.fits' % seed

def get_mog_dir():
    """path to Mixture of Gaussian directory, containing the fitted params"""
    return os.path.join(os.path.dirname(__file__),
                        '../../','etc')

def get_py_version():
    return '%s%s' % (sys.version_info[0],
                     sys.version_info[1])


def draw_points(radec,unique_ids,obj='star',seed=1,
                outdir='./'):
    """
	Args:
		radec: dict with keys ra1,ra2,dec1,dec2
			the ra,dec limits for the sample
		unique_ids: list of unique integers for each draw
		obj: star,elg,lrg,qso
		seed: to initialize random number generator
		outdir: dir to write randoms to

	Returns:
		Nothing, but write a fits_table containing the unique id, ra, dec
			and color + redshift + morphology info for each source
	"""
    print('entered draw_points')
    ndraws= len(unique_ids)
    random_state= np.random.RandomState(seed)
    ra,dec= get_radec(radec,ndraws=ndraws,random_state=random_state)
    # Load joint sample
    sample_5d_10k=fits_table(os.path.join(get_mog_dir(),
                                'elg_sample_5dim_10k.fits'))
    sample_5d_10k= fits2pandas(sample_5d_10k)
    tree = spatial.KDTree(sample_5d_10k['redshift'].values.reshape(-1,1))
    # Sample from n(z) and take the nearest z in joint sample
    model= GaussianMixtureModel.load(name=obj+'_nz',indir=get_mog_dir(),
                                     py=get_py_version(),is1D=True)
    redshifts= model.sample(ndraws)
    _,ind= tree.query(redshifts)
    boot= sample_5d_10k.iloc[ind]

    T=fits_table()
    T.set('id',unique_ids)
    # PSQL "integer" is 4 bytes
    for key in ['id']:
        T.set(key, T.get(key).astype(np.int32))
    T.set('ra',ra)
    T.set('dec',dec)
    #redshifts change each draw, not sample_5d_10k's redshift
    T.set('redshift',redshifts) 
    T.set('id_5d10k_sample',boot['id'].values)
    for col in ['g', 'r','z','rhalf']: 
        T.set(col,boot[col].values)
    # fixed priors
    if obj in ['elg','lrg']:
        T.set('n',np.ones(ndraws))
        T.set('ba', np.random.uniform(0.2,1.,size=ndraws))
        T.set('pa', np.random.uniform(0.,180.,size=ndraws))
    # Save
    fn= os.path.join(get_sample_dir(outdir,obj),get_sample_fn(seed) )
    if os.path.exists(fn):
        os.remove(fn)
        print('Overwriting %s' % fn)
    T.writeto(fn)
    print('Wrote %s' % fn)

        
def get_parser():
    parser = argparse.ArgumentParser(description='Generate a legacypipe-compatible CCDs file from a set of reduced imaging.')
    parser.add_argument('--obj', type=str, choices=['star','elg', 'lrg', 'qso'], default=None, required=True) 
    parser.add_argument('--ra1',type=float,action='store',help='bigbox',required=True)
    parser.add_argument('--ra2',type=float,action='store',help='bigbox',required=True)
    parser.add_argument('--dec1',type=float,action='store',help='bigbox',required=True)
    parser.add_argument('--dec2',type=float,action='store',help='bigbox',required=True)
    parser.add_argument('--spacing',type=float,action='store',default=10.,help='choosing N radec pionts so points have spacingxspacing arcsec spacing',required=False)
    parser.add_argument('--ndraws',type=int,action='store',help='default space by 10x10 arcsec, number of draws for all mpi tasks',required=False)
    parser.add_argument('--outdir', type=str, default='./', help='Directory to write randoms tables to')
    parser.add_argument('--nproc', type=int, default=1, help='Number of CPUs to use.')
    parser.add_argument('--seed', type=int, default=1, help='seed for nproc=1')
    parser.add_argument('--startid', type=int, default=1, help='if generating additional randoms mid-run, will want to start from a specific id')
    return parser 

if __name__ == "__main__":
    t0 = Time()
    tbegin=t0
    print('TIMING:after-imports ',datetime.datetime.now())

    parser= get_parser()
    args = parser.parse_args()
    print('TIMING:after argparse',datetime.datetime.now())

    # Before mpi, make needed dirs
    mkdir_needed( vars(args) )

    # Write calling sequence to file
    if args.nproc > 1:
        from mpi4py.MPI import COMM_WORLD as comm
        if comm.rank == 0:
            write_calling_seq( vars(args) )
    else:
        write_calling_seq( vars(args) )
         
    radec={}
    radec['ra1']=args.ra1
    radec['ra2']=args.ra2
    radec['dec1']=args.dec1
    radec['dec2']=args.dec2
    if args.ndraws is None:
        # Number that could fill a grid with 5x5 arcsec spacing
        ndraws= int( get_area(radec)/args.spacing**2 * 3600.**2 ) + 1
    else:
        ndraws= args.ndraws
    print('ndraws= %d' % ndraws)
    unique_ids= np.arange(args.startid,ndraws+args.startid)

    # Draws per mpi task
    if args.nproc > 1:
        unique_ids= np.array_split(unique_ids,comm.size)[comm.rank] 
    t0=ptime('parse-args',t0)

    if args.nproc > 1:
        seed = comm.rank
        if comm.rank == 0:
            if not os.path.exists(args.outdir):
                os.makedirs(args.outdir)
        draw_points(radec,unique_ids,obj=args.obj, seed=seed,
                    outdir=args.outdir)
    else:
        seed= args.seed
        if not os.path.exists(args.outdir):
            os.makedirs(args.outdir)
        if not os.path.exists(args.outdir):
            os.makedirs(args.outdir)
        draw_points(radec,unique_ids,obj=args.obj, seed=seed,
                    outdir=args.outdir)
        
