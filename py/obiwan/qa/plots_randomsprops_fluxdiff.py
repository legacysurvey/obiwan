if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg') # display backend
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
from scipy.stats import norm
from scipy.optimize import leastsq
from argparse import ArgumentParser

from astrometry.util.fits import fits_table, merge_tables

import obiwan.qa.plots_common as plots

plt.rcParams['axes.labelsize'] = 14
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12

#######################
parser = ArgumentParser(description='DECaLS simulations.')
parser.add_argument('--randoms_table', required=True)
parser.add_argument('--which', choices=['cosmos','eboss','desi'],required=True)
args = parser.parse_args()

dat= fits_table(args.randoms_table)

#avg_fracin= np.mean(np.array([dat.tractor_fracin_g,
#                          dat.tractor_fracin_r,
#                          dat.tractor_fracin_z]),axis=0)
#keepFracin= avg_fracin >= 0.7

fracin_thresh= 0.2
types_of_models= ['PSF','EXP','DEV','SIMP','COMP']
if args.which == 'cosmos':
    fracin_thresh= 0.1
    types_of_models= ['PSF','EXP','DEV','REX','COMP']
keepFracin= ((dat.tractor_fracin_g > fracin_thresh) &
             (dat.tractor_fracin_r > fracin_thresh) &
             (dat.tractor_fracin_z > fracin_thresh))
isRec= (dat.obiwan_mask == 1)
#isRec= (isRec) & (keepFracin)

rz= dat.psql_r - dat.psql_z
gr= dat.psql_g - dat.psql_r
mags={}
for band in 'grz':
    mags[band]= plots.flux2mag(dat.get('tractor_flux_'+band)/\
                                 dat.get('tractor_mw_transmission_'+band))

if args.which == 'eboss':
    is_elg_input= plots.eboss_ts(dat.psql_g,rz,gr,region='ngc')
    is_elg_trac= plots.eboss_ts(mags['g'],mags['r']-mags['z'],mags['g']-mags['r'],region='ngc')
#elif args.which == 'desi':
#    is_elg_input= plots.desi_ts(dat.psql_g,rz,gr)
#    is_elg_trac= plots.desi_ts(mags['g'],mags['r']-mags['z'],mags['g']-mags['r'])
print('recovered of total: %f' % (len(dat[isRec])/len(dat)))
print('recovered keepFracin of total: %f' % (len(dat[(isRec) & (keepFracin)])/len(dat)))
if args.which != 'cosmos':
    print('elg input of total: %f' % (len(dat[is_elg_input])/len(dat)))
print('keepFracin of recovered: %f' % (len(dat[(isRec) & (keepFracin)])/len(dat[isRec])))
if args.which != 'cosmos':
    print('fracin elgs >= 0.7 of recovered elgs: %f' % (len(dat[(isRec) & (keepFracin) & (is_elg_trac)])/len(dat[(isRec) & (is_elg_input)])))

##########################


def myhist(ax,data,bins=20,color='b',normed=False,lw=2,ls='solid',label=None,
           range=None, return_h=False,alpha=1):
    kw= dict(bins=bins,color=color,normed=normed,
             histtype='step',range=range,lw=lw,ls=ls,alpha=1)
    if label:
        kw.update(label=label)
    h,bins,_=ax.hist(data,**kw)
    if return_h:
        return h

def my_step(ax,bins,height,
            lw=2,color='b',ls='solid',label=None):
    """if plt.hist returns tuple (height,bins) then this reproces that plot.

    e.g. bin centers and horizontal lines at the right place...
    """
    kw= dict(color=color,lw=lw,ls=ls)
    if label:
        kw.update(label=label)
    ax.step(bins[:-1],height,where='mid',**kw)

def mytext(ax,x,y,text, ha='left',va='center',fontsize=20,rotation=0,
           color='k',dataCoords=False):
    '''adds text in x,y units of fraction axis'''
    if dataCoords:
        ax.text(x,y,text, horizontalalignment=ha,verticalalignment=va,
                fontsize=fontsize,rotation=rotation,color=color)
    else:
        ax.text(x,y,text, horizontalalignment=ha,verticalalignment=va,
                fontsize=fontsize,rotation=rotation,color=color,
                transform=ax.transAxes)

class getDepth(object):
    def __init__(self):
        self.desi= dict(g=24.0,
                        r=23.4,
                        z=22.5)
        self.eboss_ngc= dict(g=22.9,
                             r=self.desi['r'],
                             z=self.desi['z'])
        self.eboss_sgc= dict(g=22.825,
                             r=self.desi['r'],
                             z=self.desi['z'])

def grz_hist_input_noise_ext(dat,fn='grz_hist_input_noise_ext.png',
                             glim=None,rlim=None,zlim=None):
    fig,axes=plt.subplots(3,1,figsize=(5,12))
    plt.subplots_adjust(hspace=0.1,wspace=0.2)
    xlim= dict(g=glim,
               r=rlim,
               z=zlim)

    kw_hist= dict(normed=False)
    for ax,band in zip(axes,'grz'):
        flux=dict(input_noise_ext= dat.get(band+'flux'),
                  intput_noise= dat.get(band+'flux')/\
                              dat.get('mw_transmission_'+band),
                  input= plots.mag2flux(dat.get('psql_'+band))
                 )
        bins= np.linspace(xlim[band][0],xlim[band][1],num=30)
        for key,color in zip(sorted(list(flux.keys()),key=lambda x:len(x)),
                             'bgk'):
            mag= plots.flux2mag(flux[key])
            myhist(ax,mag,bins=bins,range=xlim[band],
                   color=color,label=key,**kw_hist)
        plots.mytext(ax,0.5,0.05,band,fontsize=14)
        ylab=ax.set_ylabel('Number')

    xlab=axes[-1].set_xlabel('True AB mag')
    leg=axes[0].legend(loc=(0,1.01),ncol=3)
    plt.savefig(fn,bbox_extra_artists=[xlab,ylab,leg], bbox_inches='tight')
    plt.close()
    print('Wrote %s' % fn)

def grz_hist_input_ext(dat,fn='grz_hist_input_ext.png',
                       glim=None,rlim=None,zlim=None):
    fig,axes=plt.subplots(3,1,figsize=(5,12))
    plt.subplots_adjust(hspace=0.1,wspace=0.2)
    xlim= dict(g=glim,
               r=rlim,
               z=zlim)

    kw_hist= dict(normed=False)
    for ax,band in zip(axes,'grz'):
        flux=dict(input_ext= dat.get(band+'flux'),
                  input= dat.get(band+'flux')/\
                              dat.get('mw_transmission_'+band))
        bins= np.linspace(xlim[band][0],xlim[band][1],num=30)
        for key,color in zip(sorted(list(flux.keys()),key=lambda x:len(x)),
                             'bgk'):
            mag= plots.flux2mag(flux[key])
            myhist(ax,mag,bins=bins,range=xlim[band],
                   color=color,label=key.replace('_','+'),**kw_hist)
        plots.mytext(ax,0.5,0.05,band,fontsize=14)
        ylab=ax.set_ylabel('Number')

    xlab=axes[-1].set_xlabel('True AB mag')
    leg=axes[0].legend(loc=(0,1.01),ncol=2)
    plt.savefig(fn,bbox_extra_artists=[xlab,ylab,leg], bbox_inches='tight')
    plt.close()
    print('Wrote %s' % fn)

def grz_hist_input_ext_separate_panels(dat,fn='grz_hist_input_ext.png',
                                       glim=None,rlim=None,zlim=None):
    xlim= dict(g=glim,
               r=rlim,
               z=zlim)

    kw_hist= dict(normed=False)
    for band in 'grz':
        savefn= fn.replace('.png','_%s.png' % band)
        fig,ax=plt.subplots()
        flux=dict(input_ext= dat.get(band+'flux'),
                  input= dat.get(band+'flux')/\
                              dat.get('mw_transmission_'+band))
        bins= np.linspace(xlim[band][0],xlim[band][1],num=30)
        for key,color in zip(sorted(list(flux.keys()),key=lambda x:len(x)),
                             'bgk'):
            mag= plots.flux2mag(flux[key])
            myhist(ax,mag,bins=bins,range=xlim[band],
                   color=color,label=key.replace('_','+'),**kw_hist)
        plots.mytext(ax,0.5,0.05,band,fontsize=14)
        ylab=ax.set_ylabel('Number')
        xlab=ax.set_xlabel('True AB mag')
        if band == 'g':
            leg=ax.legend(loc=(0,1.01),ncol=2)
            plt.savefig(savefn,bbox_extra_artists=[xlab,ylab,leg], bbox_inches='tight')
        else:
            plt.savefig(savefn,bbox_extra_artists=[xlab,ylab], bbox_inches='tight')
        plt.close()
        print('Wrote %s' % savefn)



def grz_hist_input_rec(dat,fn='grz_hist_input_rec.png',
                       glim=(21.5,23.25),rlim=(21.5,23.25),zlim=(19.5,22.5)):
    fig,axes=plt.subplots(3,1,figsize=(5,12))
    plt.subplots_adjust(hspace=0.2,wspace=0.2)
    xlim= dict(g=glim,
               r=rlim,
               z=zlim)

    kw_hist= dict(normed=False)
    keep= (isRec) & (keepFracin)
    for ax,band in zip(axes,'grz'):
        mag= plots.flux2mag(dat.get(band+'flux')/\
                                dat.get('mw_transmission_'+band))
        bins=np.linspace(xlim[band][0],xlim[band][1],num=30)
        myhist(ax,mag,bins=bins,
               color='b',label='input',**kw_hist)
        myhist(ax,mag[keep],bins=bins,
               color='g',label='recovered',**kw_hist)
        xlab=ax.set_xlabel('True mag %s' % band)
        ylab=ax.set_ylabel('Number')

    leg=axes[0].legend(loc=(0,1.01),ncol=2)
    plt.savefig(fn,bbox_extra_artists=[xlab,ylab,leg], bbox_inches='tight')
    plt.close()
    print('Wrote %s' % fn)

def grz_hist_by_type(dat,fn='grz_hist_by_type.png',x_ivar=0,
                     glim=(21.6,23),rlim=(20.75,22.5),zlim=(19.5,22)):

    x_ivar=0
    x_var= ['true_mag','galdepth','redshift'][x_ivar]
    types= np.char.strip(dat.get('tractor_type'))

    kw_hist=dict(bins=30,normed=True)

    figs,axes= plt.subplots(3,2,figsize=(10,9))
    plt.subplots_adjust(hspace=0.3)

    ratio_area= 1.
    xlim= dict(g=glim,
               r=rlim,
               z=zlim)
    ylim= dict(g=(0,3.7),
               r=(0,2.3),
               z=(0,1.5))
    for row,band in zip(range(3),'grz'):
        if x_var == 'true_mag':
            _x_var= plots.flux2mag(dat.get(band+'flux')/\
                                   dat.get('mw_transmission_'+band))
            xlab= '%s (true mag)' % band
        elif x_var == 'galdepth':
            flux_for_depth= 5 / np.sqrt(dat.get('tractor_galdepth_'+band))
            _x_var= plots.flux2mag(flux_for_depth/\
                                     dat.get('mw_transmission_'+band))
            xlab= 'galdepth %s' % band
        elif x_var == 'redshift':
            _x_var= dat.psql_redshift
            xlab= 'redshift'
        #isPostiveFlux= ((np.isfinite(dmag)) &
        #                (np.isfinite(true_mag)))
        #isPostiveFlux= np.ones(len(dmag),bool)
        #print('true_mag=',true_mag[isPostiveFlux],'trac_mag=',dmag[isPostiveFlux])
        for col,use_types,colors in zip(range(2),[['EXP','DEV','SIMP'],
                                                  ['SIMP','PSF','EXP']],['bgm','mcb']):
            xlabel=axes[row,col].set_xlabel(xlab)
            for typ,color in zip(use_types,colors):
                keep= (isRec) & (keepFracin) & (types == typ)
                if len(_x_var[keep]) > 0:
                    if (col == 1 and typ == 'EXP') or (col == 0 and typ == 'SIMP'):
                        #alpha= 0.25
                        lw= 0.5
                    else:
                        lw=2
                    kw_hist.update(lw=lw) #alpha=alpha)
                    myhist(axes[row,col],_x_var[keep],color=color,label=typ,range=xlim[band],**kw_hist)

        ylab='Number'
        if kw_hist['normed']:
            ylab='PDF'
        ylabel=axes[row,0].set_ylabel(ylab)
    for col in range(2):
        axes[0,col].legend(loc='upper left',ncol=1,fontsize=10)
    for row,band in zip(range(3),'grz'):
        for col in range(2):
            axes[row,col].set_xlim(xlim[band])
            axes[row,col].set_ylim(ylim[band])
    plt.savefig(fn,bbox_extra_artists=[xlabel,ylabel], bbox_inches='tight')
    plt.close()
    print('Wrote %s' % fn)

def grz_hist_elg_notelg(dat,fn='grz_hist_elg_notelg.png',
                        glim=(21.6,23),rlim=(20.75,22.5),zlim=(19.5,22)):
    fig,axes=plt.subplots(3,1,figsize=(6,12))
    plt.subplots_adjust(hspace=0.3)
    xlim= dict(g=glim,
               r=rlim,
               z=zlim)

    kw=dict(normed=True)
    for ax,band in zip(axes,'grz'):
        kw.update(bins=np.linspace(xlim[band][0],xlim[band][1],num=30),
                  range=xlim[band])
        myhist(ax,dat.get('psql_'+band)[is_elg_input],color='b',label='input, ELG',**kw)
        myhist(ax,dat.get('psql_'+band)[~is_elg_input],color='g',label='input, not ELG',**kw)

    for ax,band in zip(axes,'grz'):
        xlab=ax.set_xlabel(band)
        ylab=ax.set_ylabel('PDF')
    axes[0].legend(loc='upper left')
    plt.savefig(fn,bbox_extra_artists=[xlab,ylab], bbox_inches='tight')
    plt.close()
    print('Wrote %s' % fn)


def sum_of_noise_added(dat,fn='sum_of_noise_added.png'):
    xlim= dict(g=(21.5,23.25),
               r=(20.5,23),
               z=(19.5,22.5))

    fig,axes=plt.subplots(3,1,figsize=(5,8))
    plt.subplots_adjust(hspace=0.2)
    kw_hist= dict(color='b',bins=30,normed=False)
    ylim=2e-15
    for ax,band in zip(axes,'grz'):
        flux= plots.mag2flux(dat.get('psql_'+band))
        flux_noise= dat.get(band+'flux')/\
                        dat.get('mw_transmission_'+band)
        ax.scatter(plots.flux2mag(flux_noise),flux_noise-flux,color='b')
        ax.set_ylim(-ylim,ylim)
        xlab=ax.set_xlabel(band+' mag')
        ylab=ax.set_ylabel('dflux (input - db)')
    plt.savefig(fn,bbox_extra_artists=[xlab,ylab], bbox_inches='tight')
    plt.close()
    print('Wrote %s' % fn)

def delta_dec_vs_delta_ra(dat,fn='delta_dec_vs_delta_ra.png',
                          xlim=(-1,1),ylim=(-1,1),nbins=(30,30)):
    fig,ax= plt.subplots() #figsize=(8, 5))
    keep= (isRec) & (keepFracin)
    plots.myhist2D(ax,(dat.ra[keep] - dat.tractor_ra[keep])*3600,
                      (dat.dec[keep] - dat.tractor_dec[keep])*3600,
                   xlim=xlim,ylim=ylim,nbins=nbins)
    ax.axhline(0,c='k',ls='--')
    ax.axvline(0,c='k',ls='--')
    #plt.ylim(-1.2,1.2)
    #plt.xlim(-1.2,1.2)
    xlab=plt.xlabel(r'$\Delta \, RA$ arcsec (truth - measured)')
    ylab=plt.ylabel(r'$\Delta \, Dec$ arcsec (truth - measured)')
    plt.savefig(fn,bbox_extra_artists=[xlab,ylab], bbox_inches='tight')
    plt.close()
    print('Wrote %s' % fn)



def number_per_type_input_rec_meas(dat,fn='number_per_type_input_rec_meas.png'):
    """Horizontal Barplot

    Only count by typ injected b/c confusion matrix is the right tool if
    onsidering tractor type
    """
    use_types= ['EXP','DEV']
    injected= [len(dat[dat.n == 1]),len(dat[dat.n == 4])]
    recovered= [len(dat[(isRec) & (keepFracin) & (dat.n == 1)]),
                   len(dat[(isRec) & (keepFracin) & (dat.n == 4)])]
    df= pd.DataFrame(dict(type=use_types,
                          injected=injected,
                          recovered=recovered))
    df.set_index('type',inplace=True)

    fig,ax= plt.subplots(figsize=(8, 5))
    df.plot.barh(ax=ax)
    # Add fractions
    n_tot= np.sum(injected)
    plots.mytext(ax,0.02,0.81,'%.2f' % (recovered[1]/n_tot),fontsize=12)
    plots.mytext(ax,0.02,0.68,'%.2f' % (injected[1]/n_tot),fontsize=12)
    plots.mytext(ax,0.02,0.3,'%.2f' % (recovered[0]/n_tot),fontsize=12)
    plots.mytext(ax,0.02,0.18,'%.2f' % (injected[0]/n_tot),fontsize=12)
    xlab=ax.set_xlabel('Number')
    ylab=ax.set_ylabel('type')
    plt.savefig(fn,bbox_extra_artists=[xlab,ylab], bbox_inches='tight')
    plt.close()
    print('Wrote %s' % fn)

def confusion_matrix_by_type(dat,fn='confusion_matrix_by_type.png'):
    #use_types= ['PSF','EXP','DEV','SIMP','COMP']
    use_type= types_of_models
    trac_types= np.char.strip(dat.get('tractor_type'))
    #trac_types[pd.Series(trac_types).isin(['SIMP','REX']).values]= 'EXP'

    input_types= np.array(['EXP']*len(dat))
    keep= (isRec) & (keepFracin)
    input_types[(keep) & (dat.n == 4)]= 'DEV'
    ans_types= ['EXP','DEV']

    cm= plots.create_confusion_matrix(input_types[keep],trac_types[keep],
                                      poss_ans_types= ans_types,
                                      poss_pred_types= types_of_models)
    plt.imshow(cm, interpolation='nearest', cmap=plt.cm.Blues, vmin=0,vmax=1)
    cbar=plt.colorbar()
    plt.yticks(range(len(ans_types)), ans_types)
    plt.xticks(range(len(types_of_models)), types_of_models)
    ylab=plt.ylabel('Truth')
    xlab=plt.xlabel('Tractor')
    for row in range(len(set(input_types))):
        for col in range(len(types_of_models)):
            if np.isnan(cm[row,col]):
                plt.text(col,row,'n/a',va='center',ha='center')
            elif cm[row,col] > 0.5:
                plt.text(col,row,'%.2f' % cm[row,col],va='center',ha='center',color='yellow')
            else:
                plt.text(col,row,'%.2f' % cm[row,col],va='center',ha='center',color='black')
    plt.savefig(fn, bbox_extra_artists=[xlab,ylab], bbox_inches='tight',dpi=150)
    plt.close()
    print('Wrote %s' % fn)

def hist_true_rhalf_input(dat,fn='hist_true_rhalf_input',ylims=(0,4.2)):
    figs,ax= plt.subplots()

    bins= np.linspace(0,2,num=30)
    myhist(ax,dat.rhalf,bins=bins,color='k',
           label='Injected',normed=True)

    xlab=ax.set_xlabel(r'rhalf (true)')
    #for ax,band in zip(axes,'grz'):
    #    ax.set_xlim(ylim)
    #for ax in axes:
    ax.set_ylim(ylims)
    ylab=ax.set_ylabel('PDF')
    ax.legend(loc='upper right',fontsize=10)
    plt.savefig(fn,bbox_extra_artists=[xlab,ylab], bbox_inches='tight')
    plt.close()
    print('Wrote %s' % fn)



def hist_true_rhalf_by_type(dat,fn='hist_true_rhalf_by_type'):
    use_types= ['SIMP','EXP','DEV','PSF']
    types= np.char.strip(dat.get('tractor_type'))

    ylims= (0,7)
    figs,ax= plt.subplots(2,1,figsize=(4,6))
    plt.subplots_adjust(hspace=0.2)

    bins= np.linspace(0,2,num=30)

    keep= (isRec) & (keepFracin)
    for typ,color in zip(['SIMP','PSF'],'bc'):
        myhist(ax[0],dat.rhalf[(keep) & (types == typ)],bins=bins,color=color,
               label=typ,normed=True)
    for typ,color in zip(['EXP','DEV'],'gm'):
        myhist(ax[1],dat.rhalf[(keep) & (types == typ)],bins=bins,color=color,
               label=typ,normed=True)

    #plots.mytext(ax,0.9,0.9,typ.upper(),fontsize=14)

    #isPostiveFlux= ((np.isfinite(dmag)) &
    #                (np.isfinite(true_mag)))
    #isPostiveFlux= np.ones(len(dmag),bool)
    #print('true_mag=',true_mag[isPostiveFlux],'trac_mag=',dmag[isPostiveFlux])

    xlab=ax[-1].set_xlabel(r'rhalf (true)')
    #for ax,band in zip(axes,'grz'):
    #    ax.set_xlim(ylim)
    #for ax in axes:
    for i in range(2):
        ax[i].set_ylim(ylims)
        ylab=ax[i].set_ylabel('PDF')
        ax[i].legend(loc='upper right',fontsize=10)
    plt.savefig(fn,bbox_extra_artists=[xlab,ylab], bbox_inches='tight')
    plt.close()
    print('Wrote %s' % fn)



def fraction_recovered(dat,fn='fraction_recovered.png',
                       survey_for_depth=None,
                       glim=(20,26),rlim=(20,26),zlim=(20,26)):
    assert(survey_for_depth in ['eboss_ngc','eboss_sgc','desi'])
    fig,axes=plt.subplots(3,1,figsize=(5,12))
    plt.subplots_adjust(hspace=0.1,wspace=0.2)
    xlim= dict(g=glim,
               r=rlim,
               z=zlim)

    D= getDepth()

    kw= dict(normed=False,return_vals=True)
    keep= (isRec) & (keepFracin)
    for ax,band in zip(axes,'grz'):
        mag= plots.flux2mag(dat.get(band+'flux'))
        mag_rec= mag[keep]
        n,bins= np.histogram(mag,bins=30,range=xlim[band],normed=False)
        n_rec,_= np.histogram(mag[keep],bins=bins,range=xlim[band],normed=False)
        my_step(ax,bins,n_rec.astype(float)/n)
        ax.axhline(0.5,c='k',ls='--')
        #ax.axvline(plots.getDepth().eboss_ngc(band),c='k',ls='--')
        ax.axvline(getattr(D,survey_for_depth)[band],c='g',ls='--')
    #     ax.step(bins[:-1],n_rec/n,where='mid')
        plots.mytext(ax,0.5,0.05,band,fontsize=14)
    for ax in axes:
        ylab=ax.set_ylabel('Fraction Recovered')
        ax.set_ylim(0,1)
    xlab=axes[-1].set_xlabel('True AB mag')
    axes[0].legend(loc='upper left')
    plt.savefig(fn,bbox_extra_artists=[xlab,ylab], bbox_inches='tight')
    plt.close()
    print('Wrote %s' % fn)

def redshifts_recovered(dat,fn='redshifts_recovered.png'):
    figs,ax= plt.subplots(2,1,figsize=(4,6))
    plt.subplots_adjust(hspace=0.2)

    bins= np.linspace(0,1.5,num=30)
    # top: pdf of injected all and NGC elgs vs redshift
    myhist(ax[0],dat.psql_redshift[is_elg_input],bins=bins,color='b',
           label='NGC ELG',normed=True)
    #myhist(ax[0],dat.psql_redshift[(is_elg_input) & (isRec)],
    #       bins=bins,color='m',label='recovered',normed=True)
    myhist(ax[0],dat.psql_redshift[(is_elg_input) & (isRec) & (keepFracin) & (is_elg_trac)],
           bins=bins,color='g',label='legacypipe',normed=True)
    # bottom: fraction of recovered eboss elgs that loose to tractor measurement error
    n_elg,_= np.histogram(dat.psql_redshift[(is_elg_input)],
                                   bins=bins,normed=False)
    n_elg_legacypipe,_= np.histogram(dat.psql_redshift[(is_elg_input) & (isRec) & (keepFracin) & (is_elg_trac)],
                                   bins=bins,normed=False)
    n_notelg_legacypipe,_= np.histogram(dat.psql_redshift[(~is_elg_input) & (isRec) & (keepFracin) & (is_elg_trac)],
                                   bins=bins,normed=False)
    my_step(ax[1],bins,n_elg_legacypipe/n_elg.astype(float),
            color='g',label='legacypipe')
    my_step(ax[1],bins,n_notelg_legacypipe/n_elg.astype(float),
            color='m',label='contam by legacypipe')

    xlab=ax[-1].set_xlabel(r'redshift')
    ylab=ax[0].set_ylabel('PDF')
    ylab=ax[1].set_ylabel('Fraction')
    ax[1].set_ylim(0,1)
    for i in range(2):
        ax[i].set_xlim(bins[0],bins[-1])
        ax[i].legend(loc='upper left',fontsize=8)
    ax[0].set_yscale('log')
    plt.savefig(fn,bbox_extra_artists=[xlab,ylab], bbox_inches='tight')
    plt.close()
    print('Wrote %s' % fn)




def e1_e2(dat,fn='e1_e2.png',nbins=(120,120),
          recovered=False):
    fig,ax=plt.subplots(1,2,figsize=(8,5))
    plt.subplots_adjust(wspace=0.35)

    keep= np.ones(len(dat),bool)
    if recovered:
        keep= isRec
    plots.myhist2D(ax[0],dat.psql_ba[keep],dat.psql_pa[keep],
                   xlim=(0.1,1.1),ylim=(-20,200),nbins=nbins)
    plots.myhist2D(ax[1],dat.e1[keep],dat.e2[keep],
                   xlim=(-1,1),ylim=(-1,1),nbins=nbins)

    ax[0].set_aspect(abs((ax[0].get_xlim()[1]-ax[0].get_xlim()[0])/\
                         (ax[0].get_ylim()[1]-ax[0].get_xlim()[0])))
    ax[1].set_aspect('equal')

    for i,xlab,ylab in [(0,'ba','pa'),(1,'e1','e2')]:
        xlab=ax[i].set_xlabel(xlab)
        ylab=ax[i].set_ylabel(ylab)
    if recovered:
        fn=fn.replace('.png','_recovered.png')
    else:
        fn=fn.replace('.png','_input.png')
    plt.savefig(fn,bbox_extra_artists=[xlab,ylab], bbox_inches='tight')
    plt.close()
    print('Wrote %s' % fn)

def e1_e2_separate_panels(dat,fn='e1_e2.png',nbins=(120,120),
                          recovered=False):
    fig,ax=plt.subplots()
    keep= np.ones(len(dat),bool)
    if recovered:
        keep= isRec
    plots.myhist2D(ax,dat.e1[keep],dat.e2[keep],
                   xlim=(-1,1),ylim=(-1,1),nbins=nbins)

    ax.set_aspect('equal')

    for i,xlab,ylab in [(1,'e1','e2')]:
        xlab=ax.set_xlabel(xlab)
        ylab=ax.set_ylabel(ylab)
    if recovered:
        fn=fn.replace('.png','_recovered.png')
    else:
        fn=fn.replace('.png','_input.png')
    fn=fn.replace('.png','_1panel.png')
    plt.savefig(fn,bbox_extra_artists=[xlab,ylab], bbox_inches='tight')
    plt.close()
    print('Wrote %s' % fn)



def fraction_recovered_vs_rhalf(dat,fn='fraction_recovered_vs_rhalf.png'):
    fig,ax=plt.subplots()

    xlim=(0,2.5)
    kw=dict(normed=False,range=xlim)

    n,bins= np.histogram(dat.rhalf,bins=30,**kw)
    keep= (isRec) & (keepFracin)
    n_rec,_= np.histogram(dat.rhalf[keep],bins=bins,**kw)
    my_step(ax,bins,n_rec/n)
    ax.axhline(0.5,c='k',ls='--')

    xlab=ax.set_xlabel('rhalf (arcsec)')
    ylab=ax.set_ylabel('Fraction Recovered')
    ax.set_ylim(0.0,1.)
    plt.savefig(fn,bbox_extra_artists=[xlab,ylab], bbox_inches='tight')
    plt.close()
    print('Wrote %s' % fn)

def hist_all_quantities_fracin_cut(dat,fn='hist_all_quantities_fracin_cut.png',
                                   glim=None,rlim=None,zlim=None):
    figs,ax= plt.subplots(3,2,figsize=(8,9))
    plt.subplots_adjust(hspace=0.3,wspace=0.2)
    xlim= dict(g=glim,
               r=rlim,
               z=zlim)
    # left column is grz, right column is rhalf, redshift, number n=1,n=4
    kw_hist= dict(normed=True)
    for row,band in zip(range(3),'grz'):
        bins= np.linspace(xlim[band][0],xlim[band][1],num=30)
        true_mag= plots.flux2mag(dat.get(band+'flux')/\
                                   dat.get('mw_transmission_'+band))
        myhist(ax[row,0],true_mag[(isRec) & (keepFracin)],bins=bins,color='b',
               label='fracin >= 0.7 (good)',**kw_hist)
        myhist(ax[row,0],true_mag[(isRec) & (~keepFracin)],bins=bins,color='g',
               label='fracin < 0.7 (bad)',**kw_hist)
        xlab=ax[row,0].set_xlabel(r'%s (true mag)' % band)
        ylab=ax[row,0].set_ylabel('PDF')
    # rhalf,redshift
    bins= np.linspace(0,3,num=30)
    myhist(ax[0,1],dat.rhalf[(isRec) & (keepFracin)],bins=bins,color='b',
           label='fracin >= 0.7 (good)',**kw_hist)
    myhist(ax[0,1],dat.rhalf[(isRec) & (~keepFracin)],bins=bins,color='g',
           label='fracin < 0.7 (bad)',**kw_hist)
    bins= np.linspace(0,2,num=30)
    myhist(ax[1,1],dat.psql_redshift[(isRec) & (keepFracin)],bins=bins,color='b',
           label='fracin >= 0.7 (good)',**kw_hist)
    myhist(ax[1,1],dat.psql_redshift[(isRec) & (~keepFracin)],bins=bins,color='g',
           label='fracin < 0.7 (bad)',**kw_hist)
    for row,lab in zip([0,1],['rhalf','redshift']):
        ax[row,1].set_xlabel(lab)
        ax[row,1].set_ylabel('PDF')
    # number n=1,n=4
    frac={}
    for mod,sersic in [('exp',1),
                       ('dev',4)]:
        frac[mod+'_recovered_of_true']= len(dat[(dat.n == sersic) & (isRec)]) / len(dat[dat.n == sersic])
        print('Fraction of true %s legacypipe recovers: %.2f' %\
              (mod,frac[mod+'_recovered_of_true']))
        frac[mod+'_keepfracin_and_rec_of_true']= len(dat[(dat.n == sersic) & (isRec) & (keepFracin)]) / len(dat[(dat.n == sersic)])
        print('Fraction of recovered true %s kept after frac in:  %.2f' %\
              (mod,frac[mod+'_keepfracin_and_rec_of_true']))
    #types= np.char.strip(dat.get('tractor_type'))
    #types[pd.Series(types).isin(['SIMP','REX']).values]= 'EXP'
    #use_types= ['EXP','DEV']
    #is_true_exp= (isRec) & (dat.n == 1)
    #is_true_dev= (isRec) & (dat.n == 4)
    use_types= ['EXP','DEV']
    recovered_of_true= [frac['exp_recovered_of_true'],frac['dev_recovered_of_true']]
    keepfracin_and_rec_of_true= [frac['exp_keepfracin_and_rec_of_true'],frac['dev_keepfracin_and_rec_of_true']]
    #true_frac_rem= [len(dat[(is_true_exp) & (~keepFracin)])/len(dat[is_true_exp]),
    #                 len(dat[(is_true_dev) & (~keepFracin)])/len(dat[is_true_dev])]
    #is_trac_exp= (isRec) & (types == 'EXP')
    #is_trac_dev= (isRec) & (types == 'DEV')
    #tractor_frac_rem= [len(dat[(is_trac_exp) & (~keepFracin)])/len(dat[is_trac_exp]),
    #                    len(dat[(is_trac_dev) & (~keepFracin)])/len(dat[is_trac_dev])]
    df= pd.DataFrame(dict(type=use_types,
                          fracin=keepfracin_and_rec_of_true,
                          recovered=recovered_of_true))
    df.set_index('type',inplace=True)
    axi=ax[2,1]
    df.plot.barh(ax=axi)
    # Add fractions
    #n_tot= np.sum(injected)
    # dev
    plots.mytext(axi,0.02,0.81,'%.2f' % (recovered_of_true[1]),fontsize=12)
    plots.mytext(axi,0.02,0.68,'%.2f' % (keepfracin_and_rec_of_true[1]),fontsize=12)
    # exp
    plots.mytext(axi,0.02,0.3,'%.2f' % (recovered_of_true[0]),fontsize=12)
    plots.mytext(axi,0.02,0.18,'%.2f' % (keepfracin_and_rec_of_true[0]),fontsize=12)
    xlab=axi.set_xlabel('Fraction Remaining')
    axi.set_ylabel('')
    #for i in range(2):
    #    ax[i].set_xlim(bins[0],bins[-1])
    leg=ax[0,0].legend(loc=(0,1.01),ncol=2,fontsize=10)
    plt.savefig(fn,bbox_extra_artists=[leg,xlab,ylab], bbox_inches='tight')
    plt.close()
    print('Wrote %s' % fn)

def fracin_vs_numstddev_2dhist(dat,fn='fracin_vs_numstddev_2dhist.png',
                               delta_lims=(-6,6),nbins=(30,30)):
    figs,axes= plt.subplots(3,1,figsize=(6,10))
    plt.subplots_adjust(hspace=0.4)
    for ax,band in zip(axes,'grz'):
        x= dat.get('tractor_flux_'+band) -\
                   dat.get(band+'flux')
        x *= np.sqrt(dat.get('tractor_flux_ivar_'+band))

        y= dat.get('tractor_fracin_%s' % band)

        keep= isRec
        plots.myhist2D(ax,x[keep],y[keep],
                       xlim=delta_lims,ylim=(0,1.2),nbins=nbins)
        ax.axhline(fracin_thresh,c='r',ls='dashed')
        plots.mytext(ax,0.75,fracin_thresh+0.05,'Threshold',
                     color='r',fontsize=14)
        # label by band
        plots.mytext(ax,0.9,0.9,band,fontsize=14)
        ylab=ax.set_ylabel('fracin')
    xlab= axes[-1].set_xlabel(r'$\Delta\, Flux\,/\,\sigma$ (Tractor - Truth)')
    plt.savefig(fn,bbox_extra_artists=[xlab,ylab], bbox_inches='tight')
    plt.close()
    print('Wrote %s' % fn)




def fix_for_delta_flux(dat,fn='fix_for_delta_flux.png',
                       band='z'):
    figs,ax= plt.subplots(5,1,figsize=(5,15))
    plt.subplots_adjust(hspace=0)

    keep= (isRec) & (keepFracin)

    dflux= dat.get('tractor_flux_'+band)[keep] - dat.get(band+'flux')[keep]
    rad_aper= [0.5,0.75,1.0,1.5,2.0,3.5,5.0,7.0]
    for cnt,i_aper in zip(range(5),
                          [None,5,6,7,'avg']):
        ratio_area= 1.
        #ratio_area= (1.5*dat.rhalf[isRec] / rad_aper[i_aper])**2
        if i_aper == 'avg':
            name= 'fix: avg(-aperture_resid %.1f,%.1f)' % (rad_aper[6],rad_aper[7])
            fix= ratio_area * np.average([dat.get('tractor_apflux_resid_'+band)[keep,6],
                                          dat.get('tractor_apflux_resid_'+band)[keep,7]],
                                        axis=0)
            assert(len(fix)) == len(dat[keep])
        elif i_aper is None:
            name= 'fix: None'
            fix=0
        else:
            name= 'fix: -aperture_resid %.1f' % rad_aper[i_aper]
            fix= dat.get('tractor_apflux_resid_'+band)[keep,i_aper]* ratio_area
        y= dflux - fix
        ax[cnt].scatter(dat.get(band+'flux')[keep],y,
                        alpha=0.2,s=5,c='b',label=name)
        ax[cnt].axhline(0,c='k',ls='-',lw=2)
        ax[cnt].axhline(np.median(y),c='y',lw=2,ls='--',
                        label='Median')

    for i in range(5):
    #     ax[i].set_yscale('log')
    #     ax[i].set_ylim(1e-2,2e1)
        ax[i].set_ylim(-5,5)
        ax[i].legend(loc='upper right',markerscale=3)
    for i in range(4):
        ax[i].set_xticklabels([])
    xlab=ax[-1].set_xlabel('%s flux (nanomaggies)' % band)
    ylab=ax[-1].set_ylabel(r'$\Delta\, Flux\,/\,\sigma$ (Tractor - Truth)')
    fn=fn.replace('.png','_%s.png' % band)
    plt.savefig(fn,bbox_extra_artists=[xlab,ylab], bbox_inches='tight')
    plt.close()
    print('Wrote %s' % fn)

def delta_vs_grzmag(dat,fn='_vs_grzmag.png',
                    delta=None,delta_lims=(-6,6),typ='all',
                    nbins=(30,30),
                    glim=(17,26),rlim=(17,26),zlim=(17,26),
                    percentile_lines=True):
    assert(delta in ['num_std_dev','dmag','num_std_dev_rhalf','drhalf'])
    assert(typ in ['all','PSF','SIMP','EXP','DEV','REX'])
    fn= delta+fn.replace('.png','_bytype_%s.png' % typ)

    figs,axes= plt.subplots(3,1,figsize=(6,10))
    plt.subplots_adjust(hspace=0.4)

    xlim= dict(g=glim,
               r=rlim,
               z=zlim)
    keep= (isRec) & (keepFracin)
    for ax,band in zip(axes,'grz'):
        if delta == 'num_std_dev':
            y= dat.get('tractor_flux_'+band) -\
                       dat.get(band+'flux')
            y *= np.sqrt(dat.get('tractor_flux_ivar_'+band))
            ylabel=r'$\Delta\, Flux\,/\,\sigma$ (Tractor - Truth)'
        elif delta == 'dmag':
            # Opposite subtraction order, so < 0 mean Truth is brighter
            # Just as for delta = num_std_dev
            y= plots.flux2mag(dat.get(band+'flux')) -\
                plots.flux2mag(dat.get('tractor_flux_'+band))
            ylabel=r'$\Delta\, %s$ (Truth - Tractor)' % band
            keep= (keep) & (np.isfinite(y))
        elif delta in ['num_std_dev_rhalf','drhalf']:
            assert(typ in ['SIMP','DEV','EXP','REX'])
            if typ == 'DEV':
                eff_typ= 'dev'
            else:
                eff_typ= 'exp'
            y= dat.get('tractor_shape%s_r' % eff_typ) - dat.rhalf
            ylabel=r'$\Delta\, R_{\rm{1/2}}$ (Tractor - Truth)'
            if delta == 'num_std_dev_rhalf':
                y *= np.sqrt(dat.get('tractor_shape%s_r_ivar' % eff_typ))
                ylabel=r'$\Delta\, R_{\rm{1/2}}\,/\,\sigma$ (Tractor - Truth)'

        if typ != 'all':
            types= np.char.strip(dat.get('tractor_type'))
            keep= (keep) & (types == typ)

        true_mag= plots.flux2mag(dat.get(band+'flux')/\
                                   dat.get('mw_transmission_'+band))

        bins= np.linspace(xlim[band][0],xlim[band][1],num=30)
        plots.myhist2D(ax,true_mag[keep],y[keep],
                       xlim=xlim[band],ylim=delta_lims,nbins=nbins)

        ax.axhline(0,c='r',ls='dotted')
        if percentile_lines:
            binned= plots.bin_up(true_mag[keep],y[keep],
                                 bin_minmax=xlim[band],nbins=30)
            for perc in ['q25','q50','q75']:
                kw= dict(c='y',ls='-',lw=1)
                ax.plot(binned['binc'],binned[perc],**kw)
        ylab=ax.set_ylabel(ylabel)

    for ax,band in zip(axes,'grz'):
        xlab= ax.set_xlabel('true mag %s' % band)
        ax.legend(loc='upper left',fontsize=10,markerscale=3)
    plt.savefig(fn,bbox_extra_artists=[xlab,ylab], bbox_inches='tight')
    plt.close()
    print('Wrote %s' % fn)



def gauss_model(p,x):
    return 1/np.sqrt(2*np.pi*p[0]**2) * np.exp(-x**2/(2*p[0]**2))

def num_std_dev_gaussfit_flux(dat,fn='num_std_dev_gaussfit_flux.png',
                              delta_lims= (-6,6),typ='all',
                              keep_what_put_in='all',thresh=None,
                              sub_mean= True,cut_on_fracin=False):
    assert(typ in ['all','PSF','SIMP','EXP','DEV','REX'])
    fn= fn.replace('.png','_bytype_%s.png' % typ)
    assert(keep_what_put_in in ['all','neq1','neq4','rhalfeqpt5',
                                'neq1_notrhalf','neq4_notrhalf',
                                'fracflux','fracflux_keep_bad',
                                'fracin','fracin_keep_bad',
                                'allmask','allmask_keep_bad',
                                'fracmask','fracmask_keep_bad'])
    if keep_what_put_in != 'all':
        fn= fn.replace('.png','_keepwhatputin_%s.png' % keep_what_put_in)
    if thresh:
        fn= fn.replace('.png','_%.2f.png' % thresh)
    if not sub_mean:
        fn= fn.replace('.png','_notsubmean.png')

    figs,axes= plt.subplots(3,1,figsize=(6,10))
    plt.subplots_adjust(hspace=0.4)


    keep= isRec
    if cut_on_fracin:
        keep= (keep) & (keepFracin)
    if typ != 'all':
        types= np.char.strip(dat.get('tractor_type'))
        keep= (keep) & (types == typ)
    if keep_what_put_in != 'all':
        pad=0.05
        is_rhalf= (dat.rhalf >= 0.5-pad) & (dat.rhalf <= 0.5+pad)
        if keep_what_put_in == 'neq1':
            keep= (keep) & (dat.n == 1)
        elif keep_what_put_in == 'neq4':
            keep= (keep) & (dat.n == 4)
        elif keep_what_put_in == 'fracflux':
            fraction= np.max(np.array([dat.tractor_fracflux_g,
                                       dat.tractor_fracflux_r,
                                       dat.tractor_fracflux_z]),axis=0)
            assert(len(fraction) == len(dat))
            keep= (keep) & (fraction < thresh)
        elif keep_what_put_in in ['fracin','fracin_keep_bad']:
            #frac= np.mean(np.array([dat.tractor_fracin_g,
            #                        dat.tractor_fracin_r,
            #                        dat.tractor_fracin_z]),axis=0)
            frac= ((dat.tractor_fracin_g > thresh) &
                   (dat.tractor_fracin_r > thresh) &
                   (dat.tractor_fracin_z > thresh))
            assert(len(frac) == len(dat))
            if keep_what_put_in == 'fracin':
                #keep= (keep) & (frac > thresh) #higher fractions are good sources
                keep= (keep) & (frac)
            elif keep_what_put_in == 'fracin_keep_bad':
                #keep= (keep) & (frac <= thresh)
                keep= (keep) & (~frac)
        elif keep_what_put_in in ['fracmask','fracmask_keep_bad']:
            frac= np.mean(np.array([dat.tractor_fracmasked_g,
                                    dat.tractor_fracmasked_r,
                                    dat.tractor_fracmasked_z]),axis=0)
            assert(len(frac) == len(dat))
            if keep_what_put_in == 'fracmask':
                keep= (keep) & (frac < thresh)
            elif keep_what_put_in == 'fracmask_keep_bad':
                keep= (keep) & (frac >= thresh)
        elif keep_what_put_in in ['allmask','allmask_keep_bad']:
            good= ((dat.tractor_allmask_g == 0) &
                   (dat.tractor_allmask_r == 0) &
                   (dat.tractor_allmask_z == 0))
            if keep_what_put_in == 'allmask':
                keep= (keep) & (good)
            elif keep_what_put_in == 'allmask_keep_bad':
                keep= (keep) & (~good)
        elif keep_what_put_in == 'rhalfeqpt5':
            keep= (keep) & (is_rhalf)
        elif keep_what_put_in == 'neq1_notrhalf':
            keep= (keep) & (dat.n == 1) & (~is_rhalf)
        elif keep_what_put_in == 'neq4_notrhalf':
            keep= (keep) & (dat.n == 4) & (~is_rhalf)

    for ax,band in zip(axes,'grz'):
        data_lab= 'data'
        num_std_dev= dat.get('tractor_flux_'+band) -\
                        dat.get(band+'flux')
        num_std_dev *= np.sqrt(dat.get('tractor_flux_ivar_'+band))

        if sub_mean:
            #keep= ((num_std_dev >= num_std_lims[0]) &
            #       (num_std_dev <= num_std_lims[0])
            dflux_mean= np.mean(num_std_dev[((keep) &
                                             (num_std_dev > delta_lims[0]) &
                                             (num_std_dev < delta_lims[1]))])
            #dflux_mean= np.median(num_std_dev[isRec])
            num_std_dev -= dflux_mean
            print('%s: dflux_mean=%f' % (band,dflux_mean))
            data_lab+=' minus mean (%.2f)' % dflux_mean

        bins= np.linspace(delta_lims[0],delta_lims[1],num=30)
        h=myhist(ax,num_std_dev[keep],bins=bins,color='b',
                 label=data_lab,normed=True,
                 return_h=True)

        rv = norm()
        ax.plot(bins,rv.pdf(bins),'k--',label='Standard Norm')

        errfunc = lambda p, x, y: gauss_model(p, x) - y
        p0 = [1.] # Initial guess
        binc= (bins[:-1]+bins[1:])/2
        p1, success = leastsq(errfunc, p0[:], args=(binc, h))
        assert(success != 0)
        norm_fit= norm(scale=p1[0])
        ax.plot(bins,norm_fit.pdf(bins),'k-',label=r'Fit $\sigma=$%.2f' % p1[0])

        ax.axvline(0,c='k',ls='dotted')
        plots.mytext(ax,0.9,0.9,band,fontsize=14)
        #isPostiveFlux= ((np.isfinite(dmag)) &
        #                (np.isfinite(true_mag)))
        #isPostiveFlux= np.ones(len(dmag),bool)
        #print('true_mag=',true_mag[isPostiveFlux],'trac_mag=',dmag[isPostiveFlux])

    xlab=axes[-1].set_xlabel(r'$\Delta$flux (Tractor - True) * sqrt(ivar)')
    #for ax,band in zip(axes,'grz'):
    #    ax.set_xlim(ylim)
    for ax in axes:
        ylab=ax.set_ylabel('PDF')
        ax.legend(loc='upper left',fontsize=10,markerscale=3)
    plt.savefig(fn,bbox_extra_artists=[xlab,ylab], bbox_inches='tight')
    plt.close()
    print('Wrote %s' % fn)

def num_std_dev_gaussfit_flux_separate_panels(dat,fn='num_std_dev_gaussfit_flux_separate_plots.png',
                              delta_lims= (-6,6),
                              sub_mean= True,cut_on_fracin=False,
                              ylim=(0,0.4)):
    typ='all'
    if not sub_mean:
        fn= fn.replace('.png','_notsubmean.png')

    keep= isRec
    if cut_on_fracin:
        keep= (keep) & (keepFracin)
    for band in 'grz':
        savefn= fn.replace('.png','_%s.png' % band)
        figs,ax= plt.subplots()
        data_lab= 'data'
        num_std_dev= dat.get('tractor_flux_'+band) -\
                        dat.get(band+'flux')
        num_std_dev *= np.sqrt(dat.get('tractor_flux_ivar_'+band))

        if sub_mean:
            #keep= ((num_std_dev >= num_std_lims[0]) &
            #       (num_std_dev <= num_std_lims[0])
            dflux_mean= np.mean(num_std_dev[((keep) &
                                             (num_std_dev > delta_lims[0]) &
                                             (num_std_dev < delta_lims[1]))])
            #dflux_mean= np.median(num_std_dev[isRec])
            num_std_dev -= dflux_mean
            print('%s: dflux_mean=%f' % (band,dflux_mean))
            data_lab+=' minus mean (%.2f)' % dflux_mean

        bins= np.linspace(delta_lims[0],delta_lims[1],num=30)
        h=myhist(ax,num_std_dev[keep],bins=bins,color='b',
                 label=data_lab,normed=True,
                 return_h=True)

        rv = norm()
        ax.plot(bins,rv.pdf(bins),'k--',label='Standard Norm')

        errfunc = lambda p, x, y: gauss_model(p, x) - y
        p0 = [1.] # Initial guess
        binc= (bins[:-1]+bins[1:])/2
        p1, success = leastsq(errfunc, p0[:], args=(binc, h))
        assert(success != 0)
        norm_fit= norm(scale=p1[0])
        ax.plot(bins,norm_fit.pdf(bins),'k-',label=r'Fit $\sigma=$%.2f' % p1[0])

        ax.axvline(0,c='k',ls='dotted')
        plots.mytext(ax,0.9,0.9,band,fontsize=14)
        #isPostiveFlux= ((np.isfinite(dmag)) &
        #                (np.isfinite(true_mag)))
        #isPostiveFlux= np.ones(len(dmag),bool)
        #print('true_mag=',true_mag[isPostiveFlux],'trac_mag=',dmag[isPostiveFlux])
        ax.set_ylim(ylim)
        xlab=ax.set_xlabel(r'$\Delta$flux (Tractor - True) * sqrt(ivar)')
        ylab=ax.set_ylabel('PDF')
        ax.legend(loc='upper left',fontsize=10,markerscale=3)
        plt.savefig(savefn,bbox_extra_artists=[xlab,ylab], bbox_inches='tight')
        plt.close()
        print('Wrote %s' % savefn)


def num_std_dev_gaussfit_rhalf(dat,fn='num_std_dev_gaussfit_rhalf.png',
                               delta_lims= (-6,6),numbins=30,typ=None,
                               sub_mean= False,sub_bin_at_max=False):
    assert(typ != 'PSF') # psfsize_grz does not have ivar info
    assert(typ in ['SIMP','EXP','DEV','REX'])
    fn= fn.replace('.png','_bytype_%s.png' % typ)
    if sub_mean:
        fn= fn.replace('.png','_submean.png')

    figs,ax= plt.subplots() #figsize=(6,6))
    #plt.subplots_adjust(hspace=0.4)

    #for ax,typ in zip(axes,['exp','dev']):
    data_lab= 'data'
    #isType= types == typ.upper()
    if typ == 'DEV':
        eff_typ= 'dev'
    else:
        eff_typ= 'exp'
    rhalf= dat.get('tractor_shape%s_r' % eff_typ)
    num_std_dev= rhalf - dat.rhalf
    num_std_dev *= np.sqrt(dat.get('tractor_shape%s_r_ivar' % eff_typ))
    #keep= (np.isfinite(num_std_dev)) #num_std_dev= num_std_dev[isType]
    keep= (isRec) & (keepFracin)
    if typ != 'all':
        types= np.char.strip(dat.get('tractor_type'))
        #types[pd.Series(types).isin(['SIMP','REX']).values]= 'EXP'
        keep= (keep) & (types == typ)

    if sub_mean:
        #keep= ((num_std_dev >= num_std_lims[0]) &
        #       (num_std_dev <= num_std_lims[0])
        dflux_mean= np.mean(num_std_dev[((keep) &
                                         (num_std_dev > delta_lims[0]) &
                                         (num_std_dev < delta_lims[1]))])
        #dflux_mean= np.median(num_std_dev[isRec])
        num_std_dev -= dflux_mean
        print('%s: dflux_mean=%f' % (band,dflux_mean))
        data_lab+=' minus mean (%.2f)' % dflux_mean
    elif sub_bin_at_max:
        bins= np.linspace(delta_lims[0],delta_lims[1],num=numbins)
        h,bins=np.histogram(num_std_dev[keep],bins=bins)
        binc= (bins[:-1] + bins[1:])/2
        bin_at_max= binc[np.argmax(h)]
        num_std_dev -= bin_at_max
        data_lab+=' minus bin_at_max (%.2f)' % bin_at_max

    bins= np.linspace(delta_lims[0],delta_lims[1],num=numbins)
    h=myhist(ax,num_std_dev[keep],bins=bins,color='b',
             label=data_lab,normed=True,
             return_h=True)

    rv = norm()
    ax.plot(bins,rv.pdf(bins),'k--',label='Standard Norm')

    errfunc = lambda p, x, y: gauss_model(p, x) - y
    p0 = [1.] # Initial guess
    binc= (bins[:-1]+bins[1:])/2
    p1, success = leastsq(errfunc, p0[:], args=(binc, h))
    assert(success != 0)
    norm_fit= norm(scale=p1[0])
    ax.plot(bins,norm_fit.pdf(bins),'k-',label=r'Fit $\sigma=$%.2f' % p1[0])

    ax.axvline(0,c='k',ls='dotted')
    #plots.mytext(ax,0.9,0.9,typ.upper(),fontsize=14)

    #isPostiveFlux= ((np.isfinite(dmag)) &
    #                (np.isfinite(true_mag)))
    #isPostiveFlux= np.ones(len(dmag),bool)
    #print('true_mag=',true_mag[isPostiveFlux],'trac_mag=',dmag[isPostiveFlux])

    plots.mytext(ax,0.9,0.9,typ,fontsize=14)
    xlab=ax.set_xlabel(r'$\Delta$rhalf (Tractor - True) * sqrt(ivar)')
    #for ax,band in zip(axes,'grz'):
    #    ax.set_xlim(ylim)
    #for ax in axes:
    ylab=ax.set_ylabel('PDF')
    ax.legend(loc='upper left',fontsize=10,markerscale=3)
    plt.savefig(fn,bbox_extra_artists=[xlab,ylab], bbox_inches='tight')
    plt.close()
    print('Wrote %s' % fn)

def residual_gaussfit_rhalf(dat,fn='residual_gaussfit_rhalf.png',
                            delta_lims= (-6,6),typ=None,
                            sub_mean= False,sub_bin_at_max=False):
    """created for typ PSF b/c there is no psfsize_grz_ivar, so cannot compute num_std_dev"""
    assert(typ in ['PSF','SIMP','EXP','DEV','REX'])
    fn= fn.replace('.png','_bytype_%s.png' % typ)
    if sub_mean:
        fn= fn.replace('.png','_submean.png')

    figs,ax= plt.subplots() #figsize=(6,6))
    #plt.subplots_adjust(hspace=0.4)

    #for ax,typ in zip(axes,['exp','dev']):
    data_lab= 'data'
    #isType= types == typ.upper()
    if typ == 'PSF':
        rhalf= np.mean(np.array([dat.tractor_psfsize_g,
                                 dat.tractor_psfsize_r,
                                 dat.tractor_psfsize_z]),axis=0)/2
    elif typ == 'DEV':
        rhalf= dat.get('tractor_shape%s_r' % 'dev')
    else:
        rhalf= dat.get('tractor_shape%s_r' % 'exp')
    resid= rhalf - dat.rhalf
    #keep= (np.isfinite(num_std_dev)) #num_std_dev= num_std_dev[isType]
    keep= (isRec) & (keepFracin)
    if typ != 'all':
        types= np.char.strip(dat.get('tractor_type'))
        #types[pd.Series(types).isin(['SIMP','REX']).values]= 'EXP'
        keep= (keep) & (types == typ)

    if sub_mean:
        #keep= ((num_std_dev >= num_std_lims[0]) &
        #       (num_std_dev <= num_std_lims[0])
        dflux_mean= np.mean(resid[((keep) &
                                   (resid > delta_lims[0]) &
                                   (resid < delta_lims[1]))])
        #dflux_mean= np.median(num_std_dev[isRec])
        resid -= dflux_mean
        print('%s: dflux_mean=%f' % (band,dflux_mean))
        data_lab+=' minus mean (%.2f)' % dflux_mean
    elif sub_bin_at_max:
        bins= np.linspace(delta_lims[0],delta_lims[1],num=30)
        h,bins=np.histogram(resid[keep],bins=bins)
        binc= (bins[:-1] + bins[1:])/2
        bin_at_max= binc[np.argmax(h)]
        resid -= bin_at_max
        data_lab+=' minus bin_at_max (%.2f)' % bin_at_max

    bins= np.linspace(delta_lims[0],delta_lims[1],num=30)
    h=myhist(ax,resid[keep],bins=bins,color='b',
             label=data_lab,normed=True,
             return_h=True)

    rv = norm()
    ax.plot(bins,rv.pdf(bins),'k--',label='Standard Norm')

    errfunc = lambda p, x, y: gauss_model(p, x) - y
    p0 = [1.] # Initial guess
    binc= (bins[:-1]+bins[1:])/2
    p1, success = leastsq(errfunc, p0[:], args=(binc, h))
    assert(success != 0)
    norm_fit= norm(scale=p1[0])
    ax.plot(bins,norm_fit.pdf(bins),'k-',label=r'Fit $\sigma=$%.2f' % p1[0])

    ax.axvline(0,c='k',ls='dotted')
    #plots.mytext(ax,0.9,0.9,typ.upper(),fontsize=14)

    #isPostiveFlux= ((np.isfinite(dmag)) &
    #                (np.isfinite(true_mag)))
    #isPostiveFlux= np.ones(len(dmag),bool)
    #print('true_mag=',true_mag[isPostiveFlux],'trac_mag=',dmag[isPostiveFlux])

    xlab=ax.set_xlabel(r'$\Delta$rhalf (Tractor - True)')
    #for ax,band in zip(axes,'grz'):
    #    ax.set_xlim(ylim)
    #for ax in axes:
    ylab=ax.set_ylabel('PDF')
    ax.legend(loc='upper left',fontsize=10,markerscale=3)
    plt.savefig(fn,bbox_extra_artists=[xlab,ylab], bbox_inches='tight')
    plt.close()
    print('Wrote %s' % fn)


def num_std_dev_gaussfit_e1_e2(dat,fn='num_std_dev_gaussfit_e1_e2.png',
                               delta_lims= (-6,6),ylim=(0,0.4),typ=None,
                               sub_mean= True):
    assert(typ in ['exp','dev','simp','rex'])
    fn= fn.replace('.png','_%s.png' % typ.upper())
    types= np.char.strip(dat.get('tractor_type'))
    #types[pd.Series(types).isin(['SIMP','REX']).values]= 'EXP'

    figs,ax= plt.subplots()

    keep= (isRec) & (keepFracin) & (types == typ.upper())
    bins= np.linspace(delta_lims[0],delta_lims[1],num=30)
    rv = norm()
    ax.plot(bins,rv.pdf(bins),'k--',label='Standard Norm')
    ax.axvline(0,c='k',ls='dotted')
    for delta,color in zip(['e1','e2'],'gb'):
        data_lab= 'data'
        if typ in ['simp','rex','exp']:
            trac_e= dat.get('tractor_shapeexp_%s' % delta)
            trac_ivar= dat.get('tractor_shapeexp_%s_ivar' % delta)
        elif typ in ['dev']:
            trac_e= dat.get('tractor_shapedev_%s' % delta)
            trac_ivar= dat.get('tractor_shapedev_%s_ivar' % delta)
        if delta == 'e2':
            trac_e *= -1
        num_std_dev= trac_e -\
                        dat.get(delta)
        #print('delta=',delta)
        #print('straight diff=',num_std_dev[isType])
        num_std_dev *= np.sqrt(trac_ivar)
        #print('num_std_dev=',num_std_dev[isType])
        #print('length=',len(num_std_dev[isType]))
        #print('q25,med,q75 num_std_dev=',np.percentile(num_std_dev[isType],25),np.median(num_std_dev[isType]),np.percentile(num_std_dev[isType],75))
        good= (keep) & (np.isfinite(num_std_dev)) #num_std_dev= num_std_dev[isType]

        if sub_mean:
            #keep= ((num_std_dev >= num_std_lims[0]) &
            #       (num_std_dev <= num_std_lims[0])
            dflux_mean= np.mean(num_std_dev[((good) &
                                             (num_std_dev > delta_lims[0]) &
                                             (num_std_dev < delta_lims[1]))])
            #dflux_mean= np.median(num_std_dev[isRec])
            num_std_dev -= dflux_mean
            print('%s: dflux_mean=%f' % (band,dflux_mean))
            data_lab+=' minus mean (%.2f)' % dflux_mean

        h=myhist(ax,num_std_dev[good],bins=bins,color=color,
                 normed=True,return_h=True)

        errfunc = lambda p, x, y: gauss_model(p, x) - y
        p0 = [1.] # Initial guess
        binc= (bins[:-1]+bins[1:])/2
        p1, success = leastsq(errfunc, p0[:], args=(binc, h))
        assert(success != 0)
        norm_fit= norm(scale=p1[0])
        ax.plot(bins,norm_fit.pdf(bins),color=color,ls='-',
                label='%s - mean (%.2f),\n' % (delta,dflux_mean) + r'Fit $\sigma=$%.2f' % p1[0])

        #plots.mytext(ax,0.9,0.9,delta,fontsize=14)

        #isPostiveFlux= ((np.isfinite(dmag)) &
        #                (np.isfinite(true_mag)))
        #isPostiveFlux= np.ones(len(dmag),bool)
        #print('true_mag=',true_mag[isPostiveFlux],'trac_mag=',dmag[isPostiveFlux])

    plots.mytext(ax,0.9,0.9,typ.upper(),fontsize=14)
    xlab=ax.set_xlabel(r'$\Delta$ e1,e2 (Tractor - True) * sqrt(ivar)')
    ylab=ax.set_ylabel('PDF')
    ax.legend(loc='upper left',fontsize=10,markerscale=3)
    ax.set_ylim(ylim)
    plt.savefig(fn,bbox_extra_artists=[xlab,ylab], bbox_inches='tight')
    plt.close()
    print('Wrote %s' % fn)



def rec_lost_contam_gr_rz(dat,fn='rec_lost_contam_gr_rz.png'):
    fig,axes=plt.subplots(5,2,figsize=(10,15))
    plt.subplots_adjust(wspace=0,hspace=0)

    kw_scatter=dict(marker='.',s=20,alpha=1)
    kw_leg= dict(loc='upper left',fontsize=12,markerscale=3,frameon=False)

    good= (isRec) & (keepFracin)
    for lab,color,row,keep in [('Correct (Tractor ELG)','b',0,
                                  (good) & (is_elg_input) & (is_elg_trac)),
                               ('Contamination (Tractor ELG wrong)','g',1,
                                  (good) & (~is_elg_input) & (is_elg_trac)),
                               ('Lost (measure fails TS)','c',2,
                                  (good) & (is_elg_input) & (~is_elg_trac)),
                               ('Lost (not recovered)','m',3,
                                  (~isRec) & (is_elg_input)),
                               ('Lost (fracin)','y',4,
                                  (isRec) & (~keepFracin) & (is_elg_input))]:
        axes[row,0].scatter(dat.psql_r[keep]-dat.psql_z[keep],
                            dat.psql_g[keep]-dat.psql_r[keep],
                            c=color,label=lab,**kw_scatter)
        axes[row,1].scatter(mags['r'][keep]-mags['z'][keep],
                            mags['g'][keep]-mags['r'][keep],
                            c=color,label=lab,**kw_scatter)
        #axes[row,0].legend(**kw_leg)
        mytext(axes[row,0],0.5,0.9,lab,ha='center',fontsize=12)

    for row in range(5):
        for col in range(2):
            axes[row,col].set_xlim(0.5,2)
            axes[row,col].set_ylim(0.2,1.3)
            if row <= 2:
                axes[row,col].set_xticklabels([])
            if col == 1:
                axes[row,col].set_yticklabels([])
    for row in range(5):
        ylab=axes[row,0].set_ylabel('g-r') # (True)')
        #axes[row,1].set_ylabel('g-r (Tractor)')
    xlab=axes[-1,0].set_xlabel('r-z') # (True)')
    xlab=axes[-1,1].set_xlabel('r-z') # (Tractor)')
    title=axes[0,0].set_title('Truth',fontsize=14)
    title=axes[0,1].set_title('Tractor',fontsize=14)
    plt.savefig(fn,bbox_extra_artists=[xlab,ylab,title], bbox_inches='tight')
    plt.close()
    print('Wrote %s' % fn)

def rec_lost_contam_grz(dat,fn='rec_lost_contam_grz.png',
                        x_ivar=0):
    x_var= ['true_mag','galdepth','redshift'][x_ivar]
    kw_hist=dict(bins=30,normed=False)

    figs,axes= plt.subplots(3,1,figsize=(6,9))
    plt.subplots_adjust(hspace=0.3)

    ratio_area= 1.
    for ax,band in zip(axes,'grz'):
        if x_var == 'true_mag':
            _x_var= plots.flux2mag(dat.get(band+'flux')/\
                                   dat.get('mw_transmission_'+band))
            xlab= '%s (true mag)' % band
            xlim= dict(g=(21.6,23),
                       r=(20.75,22.5),
                       z=(19.5,22))
        elif x_var == 'galdepth':
            flux_for_depth= 5 / np.sqrt(dat.get('tractor_galdepth_'+band))
            _x_var= plots.flux2mag(flux_for_depth/\
                                     dat.get('mw_transmission_'+band))
            xlab= 'galdepth %s' % band
            xlim= dict(g=(22.5,24.5),
                       r=(22.5,24.5),
                       z=(21.5,23.5))
        elif x_var == 'redshift':
            _x_var= dat.psql_redshift
            xlab= 'redshift'
            xlim= dict(g=(0,1.5),
                       r=(0,1.5),
                       z=(0,1.5))
        #isPostiveFlux= ((np.isfinite(dmag)) &
        #                (np.isfinite(true_mag)))
        #isPostiveFlux= np.ones(len(dmag),bool)
        #print('true_mag=',true_mag[isPostiveFlux],'trac_mag=',dmag[isPostiveFlux])

        # Plot
        good= (isRec) & (keepFracin)
        for lab,color,keep in [('Correct (Tractor ELG)','b',
                                  (good) & (is_elg_input) & (is_elg_trac)),
                               ('Contamination (Tractor ELG wrong)','g',
                                  (good) & (~is_elg_input) & (is_elg_trac)),
                               ('Lost (measure fails TS)','c',
                                  (good) & (is_elg_input) & (~is_elg_trac)),
                               ('Lost (not recovered)','m',
                                  (~isRec) & (is_elg_input)),
                               ('Lost (fracin)','y',
                                  (isRec) & (~keepFracin) & (is_elg_input))]:
            myhist(ax,_x_var[keep],color=color,label=lab,range=xlim[band],**kw_hist)
        ylab='Number'
        if kw_hist['normed']:
            ylab='PDF'
        ylabel=ax.set_ylabel(ylab)
        xlabel=ax.set_xlabel(xlab)

    leg=axes[0].legend(loc=(0,1.01),ncol=2,fontsize=10,markerscale=3)
    plt.savefig(fn,bbox_extra_artists=[xlabel,ylabel,leg], bbox_inches='tight')
    plt.close()
    print('Wrote %s' % fn)


def rec_lost_contam_delta_by_type(dat,fn='rec_lost_contam_delta_by_type.png',
                                  band='g',
                                  x_ivar=0,y_ivar=0,percentile_lines=False):
    x_var= ['true_mag','galdepth','redshift'][x_ivar]
    y_var= ['dmag','chisq'][y_ivar]

    types= np.char.strip(dat.get('tractor_type'))
    use_types= ['PSF','SIMP','REX','EXP','DEV']
    if 'SIMP' in set(types):
        use_types.remove('REX')
    else:
        use_types.remove('SIMP')
    figs,axes= plt.subplots(4,1,figsize=(5,12))
    plt.subplots_adjust(hspace=0.2)

    ###
    ratio_area= 1.
    fix= ratio_area * np.average([dat.get('tractor_apflux_resid_'+band)[:,6],
                                      dat.get('tractor_apflux_resid_'+band)[:,7]],
                                    axis=0)
    assert(len(fix)) == len(dat)
    if y_var == 'dmag':
        _y_var= plots.flux2mag(dat.get('tractor_flux_'+band)) -\
                   plots.flux2mag(dat.get(band+'flux')+  fix)
        ylab= 'dmag (Tractor - True)'
        ylim=(-2,2)
    elif y_var == 'chisq':
        _y_var= dat.get('tractor_flux_'+band) -\
                   (dat.get(band+'flux')+  fix)
        _y_var *= np.sqrt(dat.get('tractor_flux_ivar_'+band))
        ylab= 'chiflux (Tractor - True) * sqrt(ivar)'
        ylim=(-10,10)

    if x_var == 'true_mag':
        _x_var= plots.flux2mag(dat.get(band+'flux')/\
                               dat.get('mw_transmission_'+band))
        xlab= 'true mag %s' % band
        xlim= dict(g=(21.6,23),
                   r=(20.75,22.5),
                   z=(19.5,22))
    elif x_var == 'galdepth':
        flux_for_depth= 5 / np.sqrt(dat.get('tractor_galdepth_'+band))
        _x_var= plots.flux2mag(flux_for_depth/\
                                 dat.get('mw_transmission_'+band))
        xlab= 'galdepth %s' % band
        xlim= dict(g=(21,24.5),
                   r=(21,24.5),
                   z=(21,24.5))
    elif x_var == 'redshift':
        _x_var= dat.psql_redshift
        xlab= 'redshift'
        xlim= dict(g=(0,1.5),
                   r=(0,1.5),
                   z=(0,1.5))
    ###

    good= (isRec) & (keepFracin)
    for ax,typ in zip(axes,use_types):
        mytext(ax,0.9,0.9,typ, fontsize=12)
        #isPostiveFlux= ((np.isfinite(dmag)) &
        #                (np.isfinite(true_mag)))
        #isPostiveFlux= np.ones(len(dmag),bool)
        #print('true_mag=',true_mag[isPostiveFlux],'trac_mag=',dmag[isPostiveFlux])

        # Plot
        xlabel=ax.set_xlabel(xlab)
        for lab,color,keep in [('lost (recovered but fail TS)','g', (good) & (is_elg_input) & (~is_elg_trac)),
                               ('Tractor ELG','b', (good) & (is_elg_input) & (is_elg_trac)),
                               ('Tractor ELG (contamiation)', 'c',(good) & (~is_elg_input) & (is_elg_trac))]:
            subset= (keep) & (types == typ)
            if percentile_lines:
                binned= plots.bin_up(_x_var[subset],_y_var[subset],
                                     bin_minmax=xlim[band],nbins=30)
                for perc in ['q25','q75']:
                    kw= dict(c=color,lw=2)
                    if perc == 'q25':
                        kw.update(label=lab)
                    ax.plot(binned['binc'],binned[perc],**kw)
            else:
                ax.scatter(_x_var[subset],_y_var[subset],
                           alpha=1,s=5,c=color,label=lab)
                #ax.scatter(true_mag[(isPostiveFlux) & (keep)],dmag[(isPostiveFlux) & (keep)],
                #           alpha=1,s=5,c=color,label=lab)

    for ax in axes:
        ax.axhline(0,c='k',ls='--')
        ax.set_ylim(ylim)
    ylabel=axes[-1].set_ylabel(ylab)
    leg=axes[0].legend(loc=(0,1.05),ncol=2,fontsize=10,markerscale=3)
    fn=fn.replace('.png','_%s.png' % band)
    plt.savefig(fn,bbox_extra_artists=[xlabel,ylabel,leg], bbox_inches='tight')
    plt.close()
    print('Wrote %s' % fn)


################
# Main
#################
if args.which == 'cosmos':
    pad=0.2
    kw_lims= dict(glim=(22-pad,24.5+pad),
                  rlim=(21.4-pad,23.9+pad),
                  zlim=(20.5-pad,23+pad))
    delta_lims=(-10,10)
    simp_or_rex='REX'
elif args.which in ['eboss','desi']:
    kw_lims= dict(glim=(21.5,23.25),
                  rlim=(20.5,23.),
                  zlim=(19.5,22.5))
    delta_lims=(-10,10)
    simp_or_rex='SIMP'

# Plots made in same order as presented in obiwan eboss paper
# Input properties
print('INPUT PROPS')
grz_hist_input_noise_ext(dat, **kw_lims)
grz_hist_input_ext(dat,**kw_lims)
grz_hist_input_ext_separate_panels(dat,**kw_lims)
e1_e2(dat,nbins=(120,120),recovered=False)
e1_e2_separate_panels(dat,nbins=(120,120),recovered=False)
sum_of_noise_added(dat)

# Identifying the num std dev peak at 1-2 sigma
# And we can remove those data points b/c they have same grz,rhalf,redshift,exp & dev frac as the other data points
print('indentifying num std dev peak at 1-2 sigma'.upper())
fracin_vs_numstddev_2dhist(dat,delta_lims=(-5,5),nbins=(30,30))
hist_all_quantities_fracin_cut(dat,**kw_lims)
num_std_dev_gaussfit_flux(dat,delta_lims= (-5,5),
                          sub_mean= False)
kw= dict(thresh=fracin_thresh,typ='all',delta_lims= (-5,5))
num_std_dev_gaussfit_flux(dat,keep_what_put_in='fracin_keep_bad',sub_mean= False,**kw)
num_std_dev_gaussfit_flux(dat,keep_what_put_in='fracin',sub_mean= True,**kw)


#########
# Rest of plots throw these data points out
print('rest of plots'.upper())
# Tractor measures input EXP much better than in put DEV
for keep_what in ['neq1','neq4']:
    num_std_dev_gaussfit_flux(dat,cut_on_fracin=True,typ='all',
                              keep_what_put_in=keep_what,
                              delta_lims= (-7,7),sub_mean= True)

delta_dec_vs_delta_ra(dat,xlim=(-1.,1.),ylim=(-1.,1.),nbins=(60,60))
grz_hist_input_rec(dat,**kw_lims)
grz_hist_by_type(dat,**kw_lims)
number_per_type_input_rec_meas(dat)
confusion_matrix_by_type(dat)
if args.which in ['eboss','desi']:
    redshifts_recovered(dat)
    fraction_recovered_vs_rhalf(dat)
fraction_recovered(dat, survey_for_depth='desi',**kw_lims)
num_std_dev_gaussfit_flux(dat,cut_on_fracin=True,typ='all',
                          delta_lims= (-5,5),sub_mean= True)
num_std_dev_gaussfit_flux_separate_panels(dat,
                          ylim=(0,0.4),delta_lims= (-5,5),
                          sub_mean= True,cut_on_fracin=True)

# rhalf measurements
hist_true_rhalf_by_type(dat)
hist_true_rhalf_input(dat)
for typ in ['PSF']:
    # Very sky distribution so sub bin at max is best
    residual_gaussfit_rhalf(dat,delta_lims= (-2,2),typ=typ,
                            sub_bin_at_max=True)
for typ in [simp_or_rex,'EXP','DEV']:
    # closer to normal, so subtract mean
    num_std_dev_gaussfit_rhalf(dat,delta_lims= delta_lims,typ=typ,
                               sub_mean=True,numbins=45)
for typ in ['exp','dev',simp_or_rex.lower()]:
    delta_vs_grzmag(dat,delta='num_std_dev_rhalf',typ=typ.upper(),
                    delta_lims=delta_lims, nbins=(60,30),**kw_lims)
    delta_vs_grzmag(dat,delta='drhalf',typ=typ.upper(),delta_lims=(-1,1),
                    nbins=(60,30),**kw_lims)
# e1,e2 measurements
for typ in ['exp','dev',simp_or_rex.lower()]:
    num_std_dev_gaussfit_e1_e2(dat,delta_lims= (-7,7),ylim=(0,0.4),
                               typ=typ,sub_mean= True)

# Combine all types
typ='all'
num_std_dev_gaussfit_flux(dat,cut_on_fracin=True,typ=typ,
                          delta_lims= (-5,5),sub_mean= True)
delta_vs_grzmag(dat,delta='num_std_dev',typ=typ,
                delta_lims=delta_lims,nbins=(60,30),**kw_lims)
delta_vs_grzmag(dat,delta='dmag',typ=typ,delta_lims=(-1,1),
                nbins=(60,30),**kw_lims)
# If want to split each of above 3 bands by type
if False:
    for typ in [simp_or_rex,'EXP','DEV','PSF']:
        num_std_dev_gaussfit_flux(dat,cut_on_fracin=True,typ=typ,
                                  delta_lims= (-5,5),sub_mean= True)
        delta_vs_grzmag(dat,delta='dmag',typ=typ,delta_lims=(-1,1),
                        nbins=(60,30),**kw_lims)
        delta_vs_grzmag(dat,delta='num_std_dev',typ=typ,delta_lims=(-10,10),
                        nbins=(60,30),**kw_lims)

# Attempt to fix mag offset by adding in sky from apflux
if False:
    for band in 'grz':
        fix_for_delta_flux(dat, band=band)

if args.which != 'cosmos':
    rec_lost_contam_gr_rz(dat)
    rec_lost_contam_grz(dat,x_ivar=0)

# FIXME: change this plot to 2D hist
for band in 'grz':
    rec_lost_contam_delta_by_type(dat,band=band,
                                  x_ivar=0,y_ivar=0,percentile_lines=False)

# Show that fracin is responsible for peak at 1-2 sigma
# And that allmask, fracmask, fracflux, rhalf ~ 0.5 are not
if False:
    for suffix in ['','_keep_bad']:
        # fracin IS RESPONSIBLE for peak at 1-2 sigma!!
        # num_std_dev for fracin < 0.7 sample and fracin >= 0.7 sample
        kw=dict(typ='all',delta_lims= (-7,7),sub_mean= False)
        for thresh in np.linspace(0.1,0.4,num=4):
            num_std_dev_gaussfit_flux(dat,keep_what_put_in='fracin'+suffix,thresh=thresh,**kw)
        # Allmask does not affect peak at 1-2 sigma
        num_std_dev_gaussfit_flux(dat,keep_what_put_in='allmask'+suffix,**kw)
        # Fracmask does not affect peak at 1-2 sigma
        for thresh in np.linspace(0.1,0.5,num=5):
            num_std_dev_gaussfit_flux(dat,keep_what_put_in='fracmask'+suffix,thresh=thresh,**kw)
        # Fracflux does not affect peak at 1-2 sigma
        for thresh in [0.01,0.5]:
            num_std_dev_gaussfit_flux(dat,keep_what_put_in='fracflux'+suffix,thresh=thresh,**kw)
    # Injected 0.45 < rhalf < 0.55 does not affect peak at 1-2 sigma
    for keep_what in ['rhalfeqpt5','neq1_notrhalf','neq4_notrhalf']:
        num_std_dev_gaussfit_flux(dat,keep_what_put_in=keep_what,**kw)
