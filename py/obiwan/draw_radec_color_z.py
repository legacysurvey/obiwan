# See LICENSE.rst for BSD 3-clause license info
# -*- coding: utf-8 -*-
"""
Uses mpi4py to draw millions of random ra,dec points and color, shape, 
redshift information from Gaussian Mixture Models. Writes a fits tables
for each mpi task.
"""

from __future__ import division, print_function

import os
import argparse
import numpy as np
import datetime
import sys
from scipy import spatial
import pandas as pd

try:
    from astrometry.util.ttime import Time
    from astrometry.util.fits import fits_table, merge_tables
except ImportError:
    pass


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
        # CANNOT use py27, the draw_points_eboss func won't work
        assert(py in ['36','35'])
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


def get_sample_fn(seed=None,startid=None):
    return 'randoms_seed_%d_startid_%d.fits' % (seed,startid)

def mog_param_dir():
    """path to Mixture of Gaussian directory, containing the fitted params"""
    return os.path.join(os.path.dirname(__file__),
                        '../../','etc')

def get_py_version():
    return '%s%s' % (sys.version_info[0],
                     sys.version_info[1])


def draw_points(radec,unique_ids,obj='star',seed=1,
                outdir='./',survey=None,startid=1):
    assert(survey in ['desi','eboss'])
    if survey == 'desi':
        return draw_points_desi(radec,unique_ids,obj=obj,seed=seed,
                                outdir=outdir,startid=startid)
    else:
        return draw_points_eboss(radec,unique_ids,obj=obj,seed=seed,
                                 outdir=outdir,startid=startid)

def draw_points_desi(radec,unique_ids,obj='star',seed=1,
                     outdir='./',startid=1):
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
    from obiwan.common import fits2pandas
    print('entered draw_points')
    ndraws= len(unique_ids)
    random_state= np.random.RandomState(seed)
    ra,dec= get_radec(radec,ndraws=ndraws,random_state=random_state)
    # Load joint sample
    sample_5d_10k=fits_table(os.path.join(outdir,
                                'elg_sample_5dim_10k.fits'))
    sample_5d_10k= fits2pandas(sample_5d_10k)
    tree = spatial.KDTree(sample_5d_10k['redshift'].values.reshape(-1,1))
    # Sample from n(z) and take the nearest z in joint sample
    model= GaussianMixtureModel.load(name=obj+'_nz',indir=mog_param_dir(),
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
    fn= os.path.join(get_sample_dir(outdir,obj),
                     get_sample_fn(seed,startid) )
    if os.path.exists(fn):
        raise IOError('fn already exists, something is wrong!, %s' % fn)
    T.writeto(fn)
    print('Wrote %s' % fn)


class EbossBox(object):
    def get_xy_pad(self,slope,pad=0):
        """Returns dx,dy"""
        theta= np.arctan(abs(slope))
        return pad*np.sin(theta), pad*np.cos(theta)
    
    def get_yint_pad(self,slope,pad=0):
        """Returns dx,dy"""
        theta= np.arctan(slope)
        return pad / np.cos(theta)

    def three_lines(self,rz,pad=0):
        slopes= np.array([-0.068,0.112, 1/(-0.555)])
        yints=  np.array([0.457,0.773,-1.901/(-0.555)])
        lines= []
        for cnt,slope,yint in zip(range(len(slopes)),slopes,yints):
            dy= 0
            #dx,dy= self.get_xy_pad(slope,pad)
            dy= self.get_yint_pad(slope,pad)
            if cnt == 0:
                dy *= -1
            #lines.append(slope*(rz-dx) + yint + dy)
            lines.append(slope*rz + yint + dy)
        return tuple(lines)
    
    def sgc_line(self,rz,pad=0):
        slope,yint= 1/0.218, -0.571/0.218
        dy=0.
        #dx,dy= self.get_xy_pad(slope,pad)
        dy= self.get_yint_pad(slope,pad)
        return slope*rz + yint + dy

    def ngc_line(self,rz,pad=0):
        slope,yint= 1/0.637, -0.399/0.637
        #dx,dy= self.get_xy_pad(slope,pad)
        dy= self.get_yint_pad(slope,pad)
        return slope*rz + yint + dy

    def SGC(self,rz, pad):
        """
        Args:
            rz: r-z
            pad: magnitudes of padding to expand TS box
        """
        d={}
        d['y1'],d['y2'],d['y3']= self.three_lines(rz,pad) 
        d['y4']= self.sgc_line(rz,pad)
        return d
    
    def NGC(self,rz, pad):
        """
        Args:
            rz: r-z
            pad: magnitudes of padding to expand TS box
        """
        d={}
        d['y1'],d['y2'],d['y3']= self.three_lines(rz,pad) 
        d['y4']= self.ngc_line(rz,pad)
        return d


def inEbossBox(rz,gr,pad=0.):
    sgc_d= EbossBox().SGC(rz,pad=pad)
    return ((gr > sgc_d['y1']) & 
            (gr < sgc_d['y2']) &
            (gr < sgc_d['y3']) &
            (gr < sgc_d['y4']))

def outside_lims_eboss(z):
    red_lims=[0.,2.]
    return ((z < red_lims[0]) | 
            (z > red_lims[1])) 


def draw_points_eboss(radec,unique_ids,obj='star',seed=1,
                      outdir='./',startid=1):
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
    # Load samples
    gmm= GaussianMixtureModel.load(name='eboss_nz_elg',indir=mog_param_dir(),
                                   py=get_py_version(),is1D=True)
    dr3dp2_exp= pd.read_csv(os.path.join(mog_param_dir(),
                            'eboss_elg_dr3deep2_EXP.csv')) # in my google drive
    dr3dp2_dev= pd.read_csv(os.path.join(mog_param_dir(),
                            'eboss_elg_dr3deep2_DEV.csv'))
    eboss_exp= pd.read_csv(os.path.join(mog_param_dir(),
                           'eboss_elg_tsspectra_EXP.csv'))
    eboss_dev= pd.read_csv(os.path.join(mog_param_dir(),
                           'eboss_elg_tsspectra_DEV.csv'))

    dr3dp2_both= pd.concat([dr3dp2_exp,dr3dp2_dev],axis='rows')
    assert(dr3dp2_both.shape[0] == dr3dp2_exp.shape[0] + dr3dp2_dev.shape[0])

    trees= dict(dr3dp2_exp= spatial.KDTree(dr3dp2_exp['redshift'].values.reshape(-1,1)),
                dr3dp2_dev= spatial.KDTree(dr3dp2_dev['redshift'].values.reshape(-1,1)),
                dr3dp2_both= spatial.KDTree(dr3dp2_both['redshift'].values.reshape(-1,1)),
                eboss_exp= spatial.KDTree(eboss_exp['redshift'].values.reshape(-1,1)),
                eboss_dev= spatial.KDTree(eboss_dev['redshift'].values.reshape(-1,1)))

    inBox=dict(dr3dp2_exp=inEbossBox(dr3dp2_exp['r'] - dr3dp2_exp['z'],
                                     dr3dp2_exp['g'] - dr3dp2_exp['r']),
               dr3dp2_dev=inEbossBox(dr3dp2_dev['r'] - dr3dp2_dev['z'],
                                     dr3dp2_dev['g'] - dr3dp2_dev['r']))

    # Draw z from n(z), if z not in [0,2] redraw
    T= fits_table()
    redshifts= gmm.sample(ndraws).reshape(-1) 
    print('redshifts=',redshifts) 

    i=0
    redraw= outside_lims_eboss(redshifts)
    num= len(redshifts[redraw])
    while num > 0:
        i+=1
        if i > 20:
            raise ValueError
        print('redrawing %d redshifts' % num)
        redshifts[redraw]= gmm.sample(num).reshape(-1)
        redraw= outside_lims_eboss(redshifts)
        num= len(redshifts[redraw])
    T.set('redshift',redshifts)

    # flip coin, 90% assign type = EXP, 10% assign type = DEV, store type
    types= np.array(['EXP']*9 + ['DEV'])
    T.set('type',types[np.random.randint(0,len(types),size=len(T))])

    # is NN redshift in the dr3_deep2 exp+dev sample, in the eboss box?
    print('len T=',len(T))
    _,i_both= trees['dr3dp2_both'].query(T.redshift.reshape(-1,1))
    print('len i_both=',len(i_both))
    #print('nn redshift=',
    inBox['dr3dp2_both']= inEbossBox(dr3dp2_both['r'].iloc[i_both] - dr3dp2_both['z'].iloc[i_both],
                                     dr3dp2_both['g'].iloc[i_both] - dr3dp2_both['r'].iloc[i_both])
    print('len inBox=',len(inBox['dr3dp2_both']))

    # Assign g,r,z,rhalf for NNs + unique id + NN redshift
    d= {}
    mag_shapes= ['g','r','z','fwhm_or_rhalf']
    for col in mag_shapes:
        d[col]= np.zeros(len(T))-1
    d['nn_redshift']= np.zeros(len(T))-1
    d['id']= np.zeros(len(T)).astype(str)
        
    # inBox, use eBOSS data
    keep= (inBox['dr3dp2_both']) & (T.type == 'EXP')
    _,i_df= trees['eboss_exp'].query(T.redshift[keep].reshape(-1,1))
    for col in mag_shapes:
        d[col][keep]= eboss_exp[col].iloc[i_df]
    d['id'][keep]= eboss_exp['sdss_id'].iloc[i_df]
    d['nn_redshift'][keep]= eboss_exp['redshift'].iloc[i_df]
        
    keep= (inBox['dr3dp2_both']) & (T.type == 'DEV')
    _,i_df= trees['eboss_dev'].query(T.redshift[keep].reshape(-1,1))
    for col in mag_shapes:
        d[col][keep]= eboss_dev[col].iloc[i_df]
    d['id'][keep]= eboss_dev['sdss_id'].iloc[i_df]
    d['nn_redshift'][keep]= eboss_dev['redshift'].iloc[i_df]

    # outBox, use DR3-Deep2 data
    keep= (~inBox['dr3dp2_both']) & (T.type == 'EXP')
    _,i_df= trees['dr3dp2_exp'].query(T.redshift[keep].reshape(-1,1))
    for col in mag_shapes:
        d[col][keep]= dr3dp2_exp[col].iloc[i_df]
    d['id'][keep]= dr3dp2_exp['tractor_id'].iloc[i_df]
    d['nn_redshift'][keep]= dr3dp2_exp['redshift'].iloc[i_df]
        
    keep= (~inBox['dr3dp2_both']) & (T.type == 'DEV')
    _,i_df= trees['dr3dp2_dev'].query(T.redshift[keep].reshape(-1,1))
    for col in mag_shapes:
        d[col][keep]= dr3dp2_dev[col].iloc[i_df]
    d['id'][keep]= dr3dp2_dev['tractor_id'].iloc[i_df]
    d['nn_redshift'][keep]= dr3dp2_dev['redshift'].iloc[i_df]

    # Add sersic n
    d['n']= np.zeros(len(T)) - 1
    d['n'][T.type == 'EXP']= 1
    d['n'][T.type == 'DEV']= 4

    for col in mag_shapes + ['nn_redshift','n']:
        if np.all(d[col] > 0) == False:
            print('FAIL: col %s has <= 0 values' % col)
            print('<= 0 values are: ',d[col][d[col] <= 0])
            raise ValueError()
    assert(np.all(pd.Series(d['id']).str.len() > 1))

    for col in mag_shapes + ['nn_redshift','n']:
        T.set(col,d[col])
    T.set('id_sample',d['id'])
    T.delete_column('type')
    T.rename('fwhm_or_rhalf','rhalf')

    # Default code for desi or eboss
    T.set('id',unique_ids)
    # PSQL "integer" is 4 bytes
    for key in ['id']:
        T.set(key, T.get(key).astype(np.int32))
    T.set('ra',ra)
    T.set('dec',dec)
    # fixed priors
    if obj in ['elg','lrg']:
        T.set('ba', np.random.uniform(0.2,1.,size=ndraws))
        T.set('pa', np.random.uniform(0.,180.,size=ndraws))
    # Save
    fn= os.path.join(get_sample_dir(outdir,obj),
                     get_sample_fn(seed,startid) )
    if os.path.exists(fn):
        raise IOError('fn already exists, something is wrong!, %s' % fn)
    T.writeto(fn)
    print('Wrote %s' % fn)

    # sanity plots
    if False:
        import seaborn as sns
        import matplotlib.pyplot as plt
        sns.distplot(T.nn_redshift - T.redshift)
        plt.savefig('delta_redshift.png')
        
        print(pd.Series(T.n).value_counts())

        cols= ['g','r','z','rhalf'] + ['redshift']
        fig,ax= plt.subplots(2,3,figsize=(12,9))
        i=-1
        for row in range(2):
            for col in range(3):
                i+=1
                if i >= len(cols):
                    continue 
                _=ax[row,col].hist(T.get(cols[i])[T.n == 1],
                                   histtype='step',normed=True,
                                   bins=30,color='b',label='EXP')
                _=ax[row,col].hist(T.get(cols[i])[T.n == 4.],
                                   histtype='step',normed=True,
                                   bins=30,color='r',label='DEV')
                ax[row,col].set_xlabel(cols[i])
        ax[1,1].legend()
        plt.savefig('hists.png')
 

        
def get_parser():
    parser = argparse.ArgumentParser(description='Generate a legacypipe-compatible CCDs file from a set of reduced imaging.')
    parser.add_argument('--survey', type=str, choices=['desi','eboss'], default=None, required=True) 
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
    parser.add_argument('--max_prev_seed', type=int, default=0, help='if generating  additional randoms need to avoid repeating a previous seed')
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
        seed = comm.rank + args.max_prev_seed
        if comm.rank == 0:
            if not os.path.exists(args.outdir):
                os.makedirs(args.outdir)
        draw_points(radec,unique_ids,obj=args.obj, seed=seed,
                    outdir=args.outdir,survey=args.survey,
                    startid=args.startid)
    else:
        seed= args.seed
        if not os.path.exists(args.outdir):
            os.makedirs(args.outdir)
        if not os.path.exists(args.outdir):
            os.makedirs(args.outdir)
        draw_points(radec,unique_ids,obj=args.obj, seed=seed,
                    outdir=args.outdir,survey=args.survey,
                    startid=args.startid)
        
