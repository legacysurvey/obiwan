import matplotlib
matplotlib.use('Agg') # display backend
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import os
import sys
from glob import glob
import numpy as np
import pandas as pd

from astrometry.util.fits import fits_table, merge_tables

import obiwan.qa.plots_common as plots

REGIONS= ['eboss_ngc','eboss_sgc','dr5','cosmos']


def footprint_outline(region):
    if region == 'eboss_ngc':
        ramin,ramax,decmin,decmax= 126,165,14,29
        x= [ramin,ramax,ramax,ramin,ramin]
        y= [decmin,decmin,decmax,decmax,decmin]
    elif region == 'eboss_sgc':
        ramin,ramax,decmin,decmax= 316,360,-2,2
        x= [ramax,ramin,ramin,ramax]
        y= [decmin,decmin,decmax,decmax]
        ramin,ramax,decmin,decmax= 0,45,-5,5
        x += [ramin,ramax,ramax,ramin,ramin]
        y += [decmax,decmax,decmin,decmin,-2]
    elif region == 'cosmos':
        lims= radec_limits(region)
        ramin,ramax= lims['ra'][0],lims['ra'][1]
        decmin,decmax= lims['dec'][0],lims['dec'][1]
        x= [ramin,ramax,ramax,ramin,ramin]
        y= [decmin,decmin,decmax,decmax,decmin]
    x,y=np.array(x), np.array(y)
    x[x >= 300] -= 360
    return x,y

def radec_limits(region):
    if region == 'eboss_ngc':
        return dict(ra=(120,170),
                    dec=(12,32))
    elif region == 'eboss_sgc':
        return dict(ra=(-50,50),
                    dec=(-6,6))
    elif region == 'cosmos':
        return dict(ra=(148.5,151.6),
                    dec=(0.4,3.4))

def cut_to_region(table,region):
    """returns a copy of table cut to region"""
    lims= radec_limits(region)
    T= table.copy()
    hiRa= T.ra >= 300
    if len(T[hiRa]) > 0:
        T.ra[hiRa] -= 360
    keep= ((T.ra >= lims['ra'][0]) & 
           (T.ra <= lims['ra'][1]) & 
           (T.dec >= lims['dec'][0]) & 
           (T.dec <= lims['dec'][1]))
    return T[keep] 



class SummaryMaps(object):
    def __init__(self,summary_fn,survey_bricks_fn):
        self.summary= fits_table(summary_fn)
        # Add ra,dec
        bricks= fits_table(survey_bricks_fn)
        keep=pd.Series(np.char.strip(bricks.brickname)).isin(np.char.strip(self.summary.brickname))
        bricks.cut(keep)
        bricks= bricks[ np.argsort(bricks.brickname)]
        self.summary= self.summary[ np.argsort(self.summary.brickname)]
        for key in ['ra','dec']:
            self.summary.set(key,bricks.get(key))

    def plot(self,region=None,subset=None):
        assert(region in REGIONS)
        keep= np.ones(len(self.summary),bool)
        if region == 'eboss_ngc':
            keep= self.summary.dec > 10
        elif region == 'eboss_sgc':
            keep= self.summary.dec < 10

        num_density= self.summary.n_injected / 0.25**2 # n/deg2
        frac_rec= self.summary.n_recovered.astype(float)/self.summary.n_injected
        d={}
        for band in 'grz':
            d['depth_'+band]= plots.flux2mag(5/np.sqrt(self.summary.get('galdepth_'+band)))

        x,y= self.summary.ra[keep],self.summary.dec[keep]
        kw=dict(x=x,y=y, figsize=(10,5),ms=15,m='o',
                region=region,subset=subset)
        if region == 'cosmos':
            kw.update(figsize=(3.6,2.9),ms=1.6e2,m='s',
                      xlim=radec_limits(region)['ra'],ylim=radec_limits(region)['dec'])
        
        kw.update(name='number_density',
                  clab='N / deg2')
        self.heatmap(num_density[keep], **kw)
        kw.update(name='fraction_recovered',
                  clab='Fraction Recovered')
        self.heatmap(frac_rec[keep],**kw)
        for band in 'grz':
            kw.update(name='galdepth_%s' % band,
                      clab='galdepth %s (not ext corr)' % band)
            newkeep= (keep) & (np.isfinite(d['depth_'+band]))
            #print('depth= ',d['depth_'+band][newkeep])
            try:
                self.heatmap(d['depth_'+band][newkeep])
            except ValueError:
                print('WARNING depth_%s FAILED, somethign about Invalid RGBA argument: 23.969765' % band)
            
    def heatmap(self,color, x=None,y=None,
                name='test',clab='Label',region=None,subset=None,
                figsize=(10,5),ms=15,m='o',
                xlim=None,ylim=None):
        fn= 'heatmap_%s_%s.png' % (name,region)
        if region == 'cosmos':
            fn= fn.replace('.png','_%s.png' % subset)
        fig,ax=plt.subplots(figsize=figsize)
        kw=dict(edgecolors='none',marker=m,s=ms,rasterized=True)
        cax= ax.scatter(x,y,c=color, **kw)
        cbar = fig.colorbar(cax) 
        
        cbar.set_label(clab)
        xlab=ax.set_xlabel('RA')
        ylab=ax.set_ylabel('Dec')
        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)
        plt.savefig(fn,bbox_extra_artists=[xlab,ylab], bbox_inches='tight')
        plt.close()
        print('Wrote %s' % fn)

class CCDsUsed(object):
    def __init__(self,ccds_used_wildcard=None):
        '''plot ra,dec of ccds used (the footprint for this run)
        
        Args:
            ccds_used_wildcard: path to ccds file or path containing an asterisk
        '''
        if "*" in ccds_used_wildcard:
            fns= glob(ccds_used_wildcard)
            assert(len(fns) > 0)
            ccds= [fits_table(fn,columns=['ra','dec','expid']) for fn in fns]
            self.ccds= merge_tables(ccds)
            del ccds
        else:
            self.ccds= fits_table(ccds_used_wildcard)

    def plot(self,savefn='ccdsused.png',
             region=None,subset=None,figsize=(10,5)):
        assert(region in REGIONS)
        savefn= savefn.replace('.png','_%s.png' % region)
        if region == 'cosmos':
            savefn= savefn.replace('.png','_%d.png' % subset)
        ccds= cut_to_region(self.ccds,region=region)
        x,y= footprint_outline(region)
        lims= radec_limits(region)
        if region == 'cosmos':
            assert(subset in [60,64,69])
            ccds= ccds[ccds.subset == subset]
        
        fig,ax= plt.subplots(figsize=figsize)
        #plt.subplots_adjust(hspace=0.2,wspace=0.2)
        # CCDs background
        if region != 'cosmos':
            ax.plot(ccds.ra,ccds.dec,'k,',lw=5,label='CCDs')
            ax.plot(x,y,'b-',lw=3)
        else:
            # ccd boundaries
            ra= dict(UL=ccds.ra0,
                     UR=ccds.ra1,
                     LR=ccds.ra2,
                     LL=ccds.ra3)
            dec= dict(UL=ccds.dec0,
                      UR=ccds.dec1,
                      LR=ccds.dec2,
                      LL=ccds.dec3)
            patches = []
            for i in range(len(ccds)):
                r= Rectangle((ra['LL'][i],dec['LL'][i]), 
                             width=ra['LR'][i]-ra['LL'][i], 
                             height=dec['UL'][i] - dec['LL'][i])
                patches.append(r)
            kw=dict(color='b',alpha=0.10)
            p = PatchCollection(patches, **kw)
            ax.add_collection(p)
            #ax.plot(x,y,'b--',lw=3)
        # Official footprint for the production run
        ax.set_xlim(lims['ra'])
        ax.set_ylim(lims['dec'])
        ax.legend(loc='lower right')
        ylab=ax.set_ylabel('Dec')
        xlab=ax.set_xlabel('RA')
        plt.savefig(savefn,bbox_extra_artists=[xlab,ylab], bbox_inches='tight')
        plt.close()
        print('Wrote %s' % savefn)

class BricksQueued(object):
    """plot ra,dec of bricks submitted to QDO queue

    Are we missing regions because didn't submit all the bricks?
    """
    def __init__(self,bricks_list_fn=None,survey_bricks_fn=None):
        '''plot ra,dec of ccds used (the footprint for this run)'''
        bricks= np.loadtxt(bricks_list_fn,dtype=str)
        self.bricks= fits_table()
        self.bricks.set('brickname',bricks)
        # Add ra,dec
        b= fits_table(survey_bricks_fn)
        keep=pd.Series(np.char.strip(b.brickname)).isin(np.char.strip(self.bricks.brickname))
        b.cut(keep)
        b= b[ np.argsort(b.brickname)]
        self.bricks= self.bricks[ np.argsort(self.bricks.brickname)]
        for key in ['ra','dec']:
            self.bricks.set(key,b.get(key))
    
    def plot(self,region=None,savefn='bricksqueued.png'):
        assert(region in REGIONS)
        savefn= savefn.replace('.png','_%s.png' % region)
        bricks= cut_to_region(self.bricks,region)
        x,y= footprint_outline(region)
        lims= radec_limits(region)

        fig,ax= plt.subplots(figsize=(10,5))
        #plt.subplots_adjust(hspace=0.2,wspace=0.2)
        # CCDs background
        ax.plot(bricks.ra,bricks.dec,'ko',markersize=1,label='Bricks')
        # Official footprint for the production run
        ax.plot(x,y,'b-',lw=3)
        ax.set_xlim(lims['ra'])
        ax.set_ylim(lims['dec'])
        ax.legend(loc='lower right')
        ylab=ax.set_ylabel('Dec')
        xlab=ax.set_xlabel('RA')
        plt.savefig(savefn,bbox_extra_artists=[xlab,ylab], bbox_inches='tight')
        plt.close()
        print('Wrote %s' % savefn)

if __name__ == "__main__":    
    from argparse import ArgumentParser
    parser = ArgumentParser(description='DECaLS simulations.')
    parser.add_argument('--region', default=None,choices=['eboss_ngc','eboss_sgc','cosmos'], required=True)
    parser.add_argument('--subset', type=int,default=None,choices=[60,64,69], required=False)
    parser.add_argument('--summary_table', default=None, required=False)
    parser.add_argument('--ccds_used_wildcard', default=None, help='path to ccds or path including wildcard for multiple ccds to merge',required=False)
    parser.add_argument('--bricks_list_fn', default=None,help='bricks.txt file', required=False)
    parser.add_argument('--survey_bricks', help='path to survey-bricks.fits.fz table', required=True)
    #parser.add_argument('--psql_radec', help='csv file from: psql -d desi -U desi_admin -h nerscdb03.nersc.gov -t -A -F"," -c "select ra,dec from obiwan_eboss_elg" > psql_radec.csv', required=True)
    args = parser.parse_args()

    if args.summary_table:
        maps= SummaryMaps(args.summary_table,args.survey_bricks)
        maps.plot(region=args.region,subset=args.subset)

    if args.ccds_used_wildcard:
        ccds= CCDsUsed(args.ccds_used_wildcard)
        kw= dict(region=args.region)
        if args.region == 'cosmos':
            kw.update(subset=args.subset)

        if args.region == 'eboss_ngc':
            kw.update(figsize=(7,5))
        elif args.region == 'eboss_sgc':
            kw.update(figsize=(10,2.5))
        elif args.region == 'cosmos':
            kw.update(figsize=(3,3))
        ccds.plot(**kw)

    if args.bricks_list_fn:
        bricks= BricksQueued(args.bricks_list_fn,args.survey_bricks)
        bricks.plot(region=args.region)
    