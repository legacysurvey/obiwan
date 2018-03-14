import matplotlib
matplotlib.use('Agg') # display backend
import matplotlib.pyplot as plt
import os
import sys
from glob import glob
import numpy as np
import pandas as pd

from astrometry.util.fits import fits_table, merge_tables

class Footprint(object):
    def __init__(self,which=None,region=None,subset=None):
        '''
        Args:
            which: eboss,dr5,cosmos
            region: ngc,sgc, subset number,None
            subset: cosmos only, integer, 60,64,69
        '''
        self.which= which
        self.region= region
        self.subset= subset

    def plot(self,ccds,psql_radec,bricks_to_run,
             uniform_randoms):
        kw=dict(ccds= self.cut_to_region(ccds),
                psql_radec= self.cut_to_region(psql_radec),
                bricks_to_run= self.cut_to_region(bricks_to_run),
                uniform_randoms= self.cut_to_region(uniform_randoms))
        #kw= dict(ccds= ccds[i['ccds']],
        #         psql_radec= psql_radec[i['psql_radec']],
        #         bricks_to_run= bricks_to_run[i['bricks_to_run']],
        #         uniform_randoms= uniform_randoms[i['uniform_randoms']])
                #ra= kw[key].ra
                #ra[split] -= 360
                #kw[key].set('ra',ra)
        
        #box= self.box_for_region()
        
        self._plot(**kw)
        #self._plot(ccds[i['ccds']],psql_radec[i['psql_radec']],
        #           bricks_to_run[i['bricks_to_run']],
        #           uniform_randoms[i['uniform_randoms']],
        #           savefn=savefn,
        #           ralim=box['ra'],declim=box['dec'])

    def _plot(self,ccds,psql_radec,bricks_to_run,
              uniform_randoms):
        fig,ax= plt.subplots(2,2,figsize=(10,10))
        plt.subplots_adjust(hspace=0.2,wspace=0.2)
        # CCDs background
        ax[0,0].plot(ccds.ra,ccds.dec,'k,',lw=5,label='CCDs')
        # Uniform randoms background
        ax[0,1].plot(uniform_randoms.ra,uniform_randoms.dec,'k,',lw=5,
                     label='Uniform randoms')
        # db randoms
        ax[1,0].plot(psql_radec.ra,psql_radec.dec,'k,',
                     label='psql randoms')
        # Bricks background
        ax[1,1].plot(bricks_to_run.ra,bricks_to_run.dec,'ko',markersize=1,label='Bricks')
        # Official footprint for the production run
        x,y= self.footprint_outline()
        lims= self.radec_limits()
        for row in range(2):
            for col in range(2):
                ax[row,col].plot(x,y,'b-',lw=3)
                ax[row,col].set_xlim(lims['ra'])
                ax[row,col].set_ylim(lims['dec'])
                ax[row,col].legend(loc='lower right')
        for i in range(2):
            ylab=ax[i,0].set_ylabel('Dec')
            xlab=ax[1,i].set_xlabel('RA')
        savefn='footprint_%s_%s.png' % (self.which,self.region)
        plt.savefig(savefn,bbox_extra_artists=[xlab,ylab], bbox_inches='tight')
        plt.close()
        print('Wrote %s' % savefn)
       

    def radec_limits(self):
        if self.which == 'eboss':
            lims= dict(ngc=dict(ra=(120,170),
                                dec=(12,32)),
                       sgc=dict(ra=(-50,50),
                                dec=(-6,6))
                      )
            return lims[self.region]
        elif self.which == 'cosmos':
            return dict(ra=(148,152),
                        dec=(0,4))

    def cut_to_region(self,table):
        lims= self.radec_limits()
        T= table.copy()
        hiRa= T.ra >= 300
        if len(T[hiRa]) > 0:
            T.ra[hiRa] -= 360
        keep= ((T.ra >= lims['ra'][0]) & 
               (T.ra <= lims['ra'][1]) & 
               (T.dec >= lims['dec'][0]) & 
               (T.dec <= lims['dec'][1]))
        if self.which == 'cosmos':
            keep= keep & (T.subset == self.subset)
        return T[keep]
         

    def footprint_outline(self):
        if self.which == 'eboss':
            if self.region == 'ngc':
                ramin,ramax,decmin,decmax= 126,165,14,29
                x= [ramin,ramax,ramax,ramin,ramin]
                y= [decmin,decmin,decmax,decmax,decmin]
            elif self.region == 'sgc':
                ramin,ramax,decmin,decmax= 316,360,-2,2
                x= [ramax,ramin,ramin,ramax]
                y= [decmin,decmin,decmax,decmax]
                ramin,ramax,decmin,decmax= 0,45,-5,5
                x += [ramin,ramax,ramax,ramin,ramin]
                y += [decmax,decmax,decmin,decmin,-2]
        elif self.which == 'cosmos':
            lims= self.radec_lims()
            ramin,ramax= lims['ra'][0],lims['ra'][1]
            decmin,decmax= lims['dec'][0],lims['dec'][1]
            x= [ramin,ramax,ramax,ramin,ramin]
            y= [decmin,decmin,decmax,decmax,decmin]
        x,y=np.array(x), np.array(y)
        x[x >= 300] -= 360
        return x,y

if __name__ == "__main__":    
    from argparse import ArgumentParser
    parser = ArgumentParser(description='DECaLS simulations.')
    parser.add_argument('--which', default=None,choices=['eboss','dr5','cosmos'], required=True)
    parser.add_argument('--path_to_ccds_used', default=None,help='./legacysurveydir_dr3', required=False)
    parser.add_argument('--ccds_used', default=None,help='single table filename', required=False)
    parser.add_argument('--psql_radec', help='csv file from: psql -d desi -U desi_admin -h nerscdb03.nersc.gov -t -A -F"," -c "select ra,dec from obiwan_eboss_elg" > psql_radec.csv', required=True)
    parser.add_argument('--bricks_to_run', help='bricks.txt file', required=True)
    parser.add_argument('--uniform_randoms', help='eboss/derived_03_05_2018/merged/randoms.fits', required=True)
    parser.add_argument('--survey_bricks', help='path to survey-bricks.fits.fz table', required=False)
    args = parser.parse_args()

    foot= Footprint(args.which)
   
    if args.ccds_used:
        ccds= fits_table(args.ccds_used)
    else:
        ccd_fns= glob(os.path.join(args.path_to_ccds_used,'survey-ccds-*-used.fits.gz'))
        assert(len(ccd_fns) > 0)
        ccds_list= [fits_table(fn,columns=['ra','dec','expid']) for fn in ccd_fns]
        ccds= merge_tables(ccds_list)
        del ccds_list

    psql_df= pd.read_csv(args.psql_radec,header=None,names=['ra','dec'])
    psql_radec=fits_table()
    for key in ['ra','dec']:
        psql_radec.set(key, psql_df[key].values)

    bricks= np.loadtxt(args.bricks_to_run,dtype=str)
    bricks_to_run= fits_table(args.survey_bricks)
    i=pd.Series(np.char.strip(bricks_to_run.brickname)).isin(bricks)
    bricks_to_run.cut(i)

    uniform_randoms= fits_table(args.uniform_randoms,
                                columns=['ra','dec'])
   
    if args.which == 'eboss':
        for region in ['ngc','sgc']:
            foot.region=region
            foot.plot(ccds,psql_radec,bricks_to_run,
                      uniform_randoms)
    elif args.which == 'cosmos':
        for subset in [60,64,69]:
            foot.subset= subset
            foot.plot(ccds,psql_radec,bricks_to_run,
                      uniform_randoms)
    

