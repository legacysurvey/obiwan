import matplotlib
matplotlib.use('Agg') # display backend
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns

from astrometry.util.fits import fits_table, merge_tables

import obiwan.qa.plots_common as plots

plt.rcParams['axes.labelsize'] = 14
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12

def eboss_ts(gmag,rz,gr,region='ngc'):
    colorCut= dict(sgc= ((gmag > 21.825) &
                         (gmag < 22.825) &
                         (-0.068 * rz + 0.457 < gr) &
                         (gr < 0.112 * rz + 0.773) &
                         (0.218 * gr + 0.571 < rz) &
                         (rz < -0.555 * gr + 1.901)),
                   ngc= ((gmag > 21.825) &
                         (gmag < 22.9) &
                         (-0.068 * rz + 0.457 < gr) &
                         (gr < 0.112 * rz + 0.773) &
                         (0.637 * gr + 0.399 < rz) &
                         (rz < -0.555 * gr + 1.901)))
    return colorCut[region]

from argparse import ArgumentParser
parser = ArgumentParser(description='DECaLS simulations.')
parser.add_argument('--randoms_table', required=True)
args = parser.parse_args()

dat= fits_table(args.randoms_table)

isRec= dat.obiwan_mask == 1
rz= dat.psql_r - dat.psql_z
gr= dat.psql_g - dat.psql_r
is_elg_input= eboss_ts(dat.psql_g,rz,gr,region='ngc')
mags={}
for band in 'grz':
    mags[band]= plots.flux2mag(dat.get('tractor_flux_'+band)/\
                                 dat.get('tractor_mw_transmission_'+band))
is_elg_trac= eboss_ts(mags['g'],mags['r']-mags['z'],mags['g']-mags['r'],region='ngc')

def myhist(ax,data,bins=20,color='b',normed=False,lw=2,ls='solid',label=None,
           range=None):
    kw= dict(bins=bins,color=color,normed=normed,
             histtype='step',range=range,lw=lw,ls=ls)
    if label:
        kw.update(label=label)
    h,bins,_=ax.hist(data,**kw)

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

def grz_hist(dat,fn='grz_hist.png'):
    fig,axes=plt.subplots(3,1,figsize=(5,12))
    plt.subplots_adjust(hspace=0.2,wspace=0.2)
    xlim= dict(g=(21.5,23.25),
               r=(20.5,23),
               z=(19.5,22.5))

    kw_hist= dict(bins=30,normed=False)
    for ax,band in zip(axes,'grz'):
        flux=dict(input_noise_ext= dat.get(band+'flux'),
                  intput_noise= dat.get(band+'flux')/\
                              dat.get('mw_transmission_'+band),
                  input= plots.mag2flux(dat.get('psql_'+band))
                 )
        for key,color in zip(sorted(list(flux.keys()),key=lambda x:len(x)),
                             'kbg'):
            mag= plots.flux2mag(flux[key])
            myhist(ax,mag,range=xlim[band],
                   color=color,label=key,**kw_hist)
        xlab=ax.set_xlabel('%s (AB mag)' % band)
        ylab=ax.set_ylabel('Number')
        
    axes[0].legend(loc='upper left')
    plt.savefig(fn,bbox_extra_artists=[xlab,ylab], bbox_inches='tight')
    plt.close()
    print('Wrote %s' % fn)

def noise_added_1(dat,fn='noise_added_1.png'):
    fig,ax=plt.subplots()

    kw_hist= dict(bins=10,normed=False)
    ylim=2e-15
    for band,color in zip('grz','gbm'):
        flux= plots.mag2flux(dat.get('psql_'+band))
        flux_noise= dat.get(band+'flux')/\
                        dat.get('mw_transmission_'+band)
        myhist(ax,flux_noise-flux,range=(-ylim,ylim),color=color,
               label=band,**kw_hist)
        ylab=ax.set_ylabel('Number')
        xlab=ax.set_xlabel('dflux (input - db)')
    ax.legend()   
    plt.savefig(fn,bbox_extra_artists=[xlab,ylab], bbox_inches='tight')
    plt.close()
    print('Wrote %s' % fn)

def noise_added_2(dat,fn='noise_added_2.png'):
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

def delta_dec_vs_delta_ra(dat,fn='delta_dec_vs_delta_ra.png'):
    plt.scatter((dat.ra[isRec] - dat.tractor_ra[isRec])*3600,
                (dat.dec[isRec] - dat.tractor_dec[isRec])*3600,
                alpha=0.2,s=5,c='b')
    plt.axhline(0,c='k',ls='--')
    plt.axvline(0,c='k',ls='--')
    plt.ylim(-1.2,1.2)
    plt.xlim(-1.2,1.2)
    xlab=plt.xlabel(r'$\Delta \, RA$ (truth - measured)')
    ylab=plt.ylabel(r'$\Delta \, Dec$ (truth - measured)')
    plt.savefig(fn,bbox_extra_artists=[xlab,ylab], bbox_inches='tight')
    plt.close()
    print('Wrote %s' % fn)
 
def fraction_recovered(dat,fn='fraction_recovered.png',
                       eboss_or_desi='eboss'):
    fig,axes=plt.subplots(3,1,figsize=(5,12))
    plt.subplots_adjust(hspace=0.2,wspace=0.2)
    xlim= dict(g=(21.5,23.25),
               r=(20.5,23),
               z=(19.5,22.5))

    D= getDepth()

    kw= dict(normed=False,return_vals=True)
    for ax,band in zip(axes,'grz'):
        mag= plots.flux2mag(dat.get(band+'flux'))
        mag_rec= mag[isRec]
        n,bins= np.histogram(mag,bins=30,range=xlim[band],normed=False)
        n_rec,_= np.histogram(mag[isRec],bins=bins,range=xlim[band],normed=False)
        my_step(ax,bins,n_rec.astype(float)/n)
        ax.axhline(0.5,c='k',ls='--')
        #ax.axvline(plots.getDepth().eboss_ngc(band),c='k',ls='--')
        ax.axvline(D.eboss_ngc[band],c='k',ls='--')
    #     ax.step(bins[:-1],n_rec/n,where='mid')
        xlab=ax.set_xlabel('%s (AB mag)' % band)
    for ax in axes:
        ylab=ax.set_ylabel('Fraction Recovered')
        ax.set_ylim(0,1)
    axes[0].legend(loc='upper left')
    plt.savefig(fn,bbox_extra_artists=[xlab,ylab], bbox_inches='tight')
    plt.close()
    print('Wrote %s' % fn)

def e1_e2_input(dat,fn='e1_e2_input.png'):
    fig,ax=plt.subplots(1,2,figsize=(8,5))
    plt.subplots_adjust(wspace=0.3)

    kw=dict(color='b',m='o',s=10.,alpha=0.25)
    plots.myscatter_open(ax[0],dat.e1,dat.e2,**kw)
    plots.myscatter_open(ax[1],dat.psql_ba,dat.psql_pa,**kw)

    ax[0].set_aspect('equal')
    ax[1].set_aspect(abs((ax[1].get_xlim()[1]-ax[1].get_xlim()[0])/\
                         (ax[1].get_ylim()[1]-ax[1].get_xlim()[0])))

    for i,xlab,ylab in [(0,'e1','e2'),(1,'ba','pa')]:
        xlab=ax[i].set_xlabel(xlab)
        ylab=ax[i].set_ylabel(ylab)
    plt.savefig(fn,bbox_extra_artists=[xlab,ylab], bbox_inches='tight')
    plt.close()
    print('Wrote %s' % fn)

 
def e1_e2_recovered(dat,fn='e1_e2_recovered.png'):
    fig,ax=plt.subplots(1,2,figsize=(8,5))
    plt.subplots_adjust(wspace=0.3)

    kw=dict(color='b',m='o',s=10.,alpha=0.25)
    plots.myscatter_open(ax[0],dat.e1[isRec],dat.e2[isRec],**kw)
    plots.myscatter_open(ax[1],dat.psql_ba[isRec],dat.psql_pa[isRec],**kw)

    ax[0].set_aspect('equal')
    ax[1].set_aspect(abs((ax[1].get_xlim()[1]-ax[1].get_xlim()[0])/\
                         (ax[1].get_ylim()[1]-ax[1].get_xlim()[0])))

    for i,xlab,ylab in [(0,'e1','e2'),(1,'ba','pa')]:
        xlabel=ax[i].set_xlabel(xlab)
        ylabel=ax[i].set_ylabel(ylab)
    plt.savefig(fn,bbox_extra_artists=[xlabel,ylabel], bbox_inches='tight')
    plt.close()
    print('Wrote %s' % fn)

def fraction_recovered_vs_rhalf(dat,fn='fraction_recovered_vs_rhalf.png'):
    fig,ax=plt.subplots()
        
    xlim=(0,2.5)
    kw=dict(normed=False,range=xlim)

    n,bins= np.histogram(dat.rhalf,bins=30,**kw)
    n_rec,_= np.histogram(dat.rhalf[isRec],bins=bins,**kw)
    my_step(ax,bins,n_rec/n)
    ax.axhline(0.5,c='k',ls='--')

    xlab=ax.set_xlabel('rhalf (arcsec)')
    ylab=ax.set_ylabel('Fraction Recovered')
    ax.set_ylim(0.0,1.)
    plt.savefig(fn,bbox_extra_artists=[xlab,ylab], bbox_inches='tight')
    plt.close()
    print('Wrote %s' % fn)

def rec_lost_contam_gr_rz_g(dat,fn='rec_lost_contam_gr_rz_g.png'):
    fig,axes=plt.subplots(3,2,figsize=(12,9))
    plt.subplots_adjust(wspace=0.2,hspace=0.4)

    kw_scatter=dict(marker='.',s=20,alpha=1)
    kw_hist=dict(range=(21.5,23.5),normed=False)

    # Row 0
    row=0
    lab,keep= 'true ELG',isRec & is_elg_input
    axes[row,0].scatter(dat.psql_r[keep]-dat.psql_z[keep],
                  dat.psql_g[keep]-dat.psql_r[keep],
                  c='k',label=lab,**kw_scatter)
    myhist(axes[row,1],dat.psql_g[keep],bins=30,color='k',label=lab,**kw_hist)

    # Rows 1 & 2
    mags={}
    for band in 'grz':
        mags[band]= plots.flux2mag(dat.get('tractor_flux_'+band)/\
                                     dat.get('tractor_mw_transmission_'+band))
    is_elg_trac= eboss_ts(mags['g'],mags['r']-mags['z'],mags['g']-mags['r'],region='ngc')

    lab,keep= 'lost (recovered but fail TS)', (isRec) & (is_elg_input) & (~is_elg_trac)
    row=1
    axes[row,0].scatter(dat.psql_r[keep]-dat.psql_z[keep],
                        dat.psql_g[keep]-dat.psql_r[keep],
                        c='g',label=lab,**kw_scatter)
    myhist(axes[row,1],dat.psql_g[keep],bins=30,color='g',label=lab,**kw_hist)
    row=2
    axes[row,0].scatter(mags['r'][keep]-mags['z'][keep],
                        mags['g'][keep]-mags['r'][keep],
                        c='g',label=lab,**kw_scatter)
    myhist(axes[row,1],mags['g'][keep],bins=30,color='g',label=lab,**kw_hist)

    lab,keep= 'Tractor ELG (correct)', (isRec) & (is_elg_input) & (is_elg_trac)
    row=1
    axes[row,0].scatter(dat.psql_r[keep]-dat.psql_z[keep],
                        dat.psql_g[keep]-dat.psql_r[keep],
                        c='b',label=lab,**kw_scatter)
    myhist(axes[row,1],dat.psql_g[keep],bins=30,color='b',label=lab,**kw_hist)
    row=2
    axes[row,0].scatter(mags['r'][keep]-mags['z'][keep],
                        mags['g'][keep]-mags['r'][keep],
                        c='b',label=lab,**kw_scatter)
    myhist(axes[row,1],mags['g'][keep],bins=30,color='b',label=lab,**kw_hist)

    lab,keep= 'Tractor ELG (contamiation)', (isRec) & (~is_elg_input) & (is_elg_trac)
    row=1
    axes[row,0].scatter(dat.psql_r[keep]-dat.psql_z[keep],
                        dat.psql_g[keep]-dat.psql_r[keep],
                        c='b',label=lab,**kw_scatter)
    myhist(axes[row,1],dat.psql_g[keep],bins=30,color='c',label=lab,**kw_hist)
    row=2
    axes[row,0].scatter(mags['r'][keep]-mags['z'][keep],
                        mags['g'][keep]-mags['r'][keep],
                        c='b',label=lab,**kw_scatter)
    myhist(axes[row,1],mags['g'][keep],bins=30,color='c',label=lab,**kw_hist)


    lab,keep= 'lost (not recovered)', (~isRec) & (is_elg_input)
    row=1
    axes[row,0].scatter(dat.psql_r[keep]-dat.psql_z[keep],
                        dat.psql_g[keep]-dat.psql_r[keep],
                        c='m',label=lab,**kw_scatter)
    myhist(axes[row,1],dat.psql_g[keep],bins=30,color='m',label=lab,**kw_hist)
    row=2
    axes[row,0].scatter(mags['r'][keep]-mags['z'][keep],
                        mags['g'][keep]-mags['r'][keep],
                        c='m',label=lab,**kw_scatter)
    myhist(axes[row,1],mags['g'][keep],bins=30,color='m',label=lab,**kw_hist)


    col=0
    for row in range(3):
        axes[row,col].set_xlim(0.5,2)
        axes[row,col].set_ylim(0.2,1.3)
    axes[0,col].set_xlabel('True r-z')
    axes[0,col].set_ylabel('True g-r')
    axes[1,col].set_xlabel('True r-z')
    axes[1,col].set_ylabel('True g-r')
    axes[2,col].set_xlabel('Tractor r-z')
    ylab=axes[2,col].set_ylabel('Tractor g-r')

    col=1
    axes[0,col].set_xlabel('True g')
    axes[1,col].set_xlabel('True g')
    xlab=axes[2,col].set_xlabel('Tractor g')
    kw_leg= dict(loc='upper left',fontsize=8)
    axes[0,col].legend(**kw_leg)
    axes[1,col].legend(**kw_leg)
    axes[2,col].legend(**kw_leg) 
    plt.savefig(fn,bbox_extra_artists=[xlab,ylab], bbox_inches='tight')
    plt.close()
    print('Wrote %s' % fn)

def rec_lost_contam_grz(dat,fn='rec_lost_contam_grz.png',
                        x_ivar=0):
    x_var= ['true_mag','galdepth','redshift'][x_ivar]
    kw_hist=dict(bins=30,normed=False)

    figs,axes= plt.subplots(3,1,figsize=(5,10))
    plt.subplots_adjust(hspace=0.4)

    ratio_area= 1. 
    for ax,band in zip(axes,'grz'):
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
        for lab,color,keep in [('lost (recovered but fail TS)','g', (isRec) & (is_elg_input) & (~is_elg_trac)),
                               ('Tractor ELG','b', (isRec) & (is_elg_input) & (is_elg_trac)),
                               ('Tractor ELG (contamiation)', 'c',(isRec) & (~is_elg_input) & (is_elg_trac))]:
            myhist(ax,_x_var[keep],color=color,label=lab,range=xlim[band],**kw_hist)
        ylab='Number'
        if kw_hist['normed']:
            ylab='PDF'
        ylabel=ax.set_ylabel(ylab)
        xlabel=ax.set_xlabel(xlab)

    leg=axes[0].legend(loc=(0,1.01),ncol=2,fontsize=12,markerscale=3)
    plt.savefig(fn,bbox_extra_artists=[xlabel,ylabel,leg], bbox_inches='tight')
    plt.close()
    print('Wrote %s' % fn)

def rec_lost_contam_by_type(dat,fn='rec_lost_contam_by_type',
                            band='g',x_ivar=0):
    x_ivar=0
    x_var= ['true_mag','galdepth','redshift'][x_ivar]
    types= np.char.strip(dat.get('tractor_type'))
    use_types= ['PSF','SIMP','REX','EXP','DEV','COMP']

    kw_hist=dict(bins=30,normed=False)

    figs,axes= plt.subplots(1,3,figsize=(12,4))
    plt.subplots_adjust(wspace=0.2)

    ratio_area= 1. 
    for ax,which in zip(axes,
                        ['lost (recovered but fail TS)','Tractor ELG',
                         'Tractor ELG (contamiation)']):
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
        if 'lost' in which:
            keep= (isRec) & (is_elg_input) & (~is_elg_trac)
        elif 'contam' in which:
            keep= (isRec) & (~is_elg_input) & (is_elg_trac)
        else:
            keep= (isRec) & (is_elg_input) & (is_elg_trac)
        ax.set_title(which)
        
        # Plot
        for typ,color in zip(use_types,'kbbgmc'):
            subset= (keep) & (types == typ)
            if len(_x_var[subset]) > 0:
                myhist(ax,_x_var[subset],color=color,label=typ,range=xlim[band],**kw_hist)
        ylab='Number'
        if kw_hist['normed']:
            ylab='PDF'
        ylabel=ax.set_ylabel(ylab)
        xlabel=ax.set_xlabel(xlab)
    axes[0].legend(loc='upper left',ncol=1,fontsize=10)
    plt.savefig(fn,bbox_extra_artists=[xlabel,ylabel], bbox_inches='tight')
    plt.close()
    print('Wrote %s' % fn)


def rec_lost_contam_delta(dat,fn='rec_lost_contam_delta.png',
                          x_ivar=0,y_ivar=0,
                          percentile_lines=False):
    x_var= ['true_mag','galdepth','redshift'][x_ivar]
    y_var= ['dmag','chisq'][y_ivar]

    figs,axes= plt.subplots(3,1,figsize=(6,10))
    plt.subplots_adjust(hspace=0.4)

    ratio_area= 1. 
    for ax,band in zip(axes,'grz'):
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
        
        #isPostiveFlux= ((np.isfinite(dmag)) &
        #                (np.isfinite(true_mag)))
        #isPostiveFlux= np.ones(len(dmag),bool)
        #print('true_mag=',true_mag[isPostiveFlux],'trac_mag=',dmag[isPostiveFlux])
        
        # Plot
        xlabel=ax.set_xlabel(xlab)
        for lab,color,keep in [('lost (recovered but fail TS)','g', (isRec) & (is_elg_input) & (~is_elg_trac)),
                               ('Tractor ELG','b', (isRec) & (is_elg_input) & (is_elg_trac)),
                               ('Tractor ELG (contamiation)', 'c',(isRec) & (~is_elg_input) & (is_elg_trac))]:
            if percentile_lines:
                binned= plots.bin_up(_x_var[keep],_y_var[keep], 
                                     bin_minmax=xlim[band],nbins=30)
                for perc in ['q25','q75']:
                    kw= dict(c=color,lw=2)
                    if perc == 'q25':
                        kw.update(label=lab)
                    ax.plot(binned['binc'],binned[perc],**kw)
            else:
                ax.scatter(_x_var[keep],_y_var[keep],
                           alpha=1,s=5,c=color,label=lab)
                #ax.scatter(true_mag[(isPostiveFlux) & (keep)],dmag[(isPostiveFlux) & (keep)],
                #           alpha=1,s=5,c=color,label=lab)
        
    for ax,band in zip(axes,'grz'):
        ax.axhline(0,c='k',ls='--')
        ax.set_ylim(ylim)
    ylabel=axes[1].set_ylabel(ylab)
    leg=axes[0].legend(loc=(0,1.01),ncol=2,fontsize=12,markerscale=3)
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

    for ax,typ in zip(axes,use_types):
        mytext(ax,0.9,0.9,typ, fontsize=12)
        #isPostiveFlux= ((np.isfinite(dmag)) &
        #                (np.isfinite(true_mag)))
        #isPostiveFlux= np.ones(len(dmag),bool)
        #print('true_mag=',true_mag[isPostiveFlux],'trac_mag=',dmag[isPostiveFlux])
        
        # Plot
        xlabel=ax.set_xlabel(xlab)
        for lab,color,keep in [('lost (recovered but fail TS)','g', (isRec) & (is_elg_input) & (~is_elg_trac)),
                               ('Tractor ELG','b', (isRec) & (is_elg_input) & (is_elg_trac)),
                               ('Tractor ELG (contamiation)', 'c',(isRec) & (~is_elg_input) & (is_elg_trac))]:
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
    plt.savefig(fn,bbox_extra_artists=[xlabel,ylabel,leg], bbox_inches='tight')
    plt.close()
    print('Wrote %s' % fn)

def rec_lost_contam_input_elg_notelg(dat,fn='rec_lost_contam_input_elg_notelg.png'):
    fig,axes=plt.subplots(3,1,figsize=(6,12))
    plt.subplots_adjust(hspace=0.3)
    xlim= dict(g=(21.6,23),
               r=(20.75,22.5),
               z=(19.5,22))

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


def rec_lost_contam_fraction(dat,fn='rec_lost_contam_fraction.png'):
    fig,axes=plt.subplots(1,2,figsize=(10,5))
    plt.subplots_adjust(wspace=0.2)
    # xlim=(0,2.5)
    xlim=None
    kw=dict(range=(0,2),normed=False)

    # Left panel
    ax= axes[0]
    n_tot,bins= np.histogram(dat.psql_redshift,bins=30,**kw)

    lab,keep= 'is ELG', (is_elg_input)
    n_elg,bins= np.histogram(dat.psql_redshift[keep],bins=bins,**kw)
    my_step(ax,bins,n_elg/n_tot,color='b',label=lab)

    lab,keep= 'is ELG, recovered', (isRec) & (is_elg_input)
    n_elg_rec,_= np.histogram(dat.psql_redshift[keep],bins=bins,**kw)
    my_step(ax,bins,n_elg_rec/n_tot,color='g',label=lab)

    # Right panel
    ax= axes[1]
    lab,keep= 'is ELG, recovered ELG', (isRec) & (is_elg_input) & (is_elg_trac)
    n_corr,_= np.histogram(dat.psql_redshift[keep],bins=bins,**kw)
    my_step(ax,bins,n_corr/n_elg_rec,color='b',label=lab)
    print('n_corr=',n_corr)

    lab,keep= 'is ELG, recovered !ELG', (isRec) & (is_elg_input) & (~is_elg_trac)
    n_lost,_= np.histogram(dat.psql_redshift[keep],bins=bins,**kw)
    my_step(ax,bins,n_lost/n_elg_rec,color='g',label=lab)
    print('n_lost=',n_lost)

    lab,keep= 'not ELG, recovered ELG', (isRec) & (~is_elg_input) & (is_elg_trac)
    n_contam,_= np.histogram(dat.psql_redshift[keep],bins=bins,**kw)
    my_step(ax,bins,n_contam/n_elg_rec,color='m',label=lab)
    print('n_contam=',n_contam)

    for ax in axes:
        ax.axhline(0.5,c='k',ls='--')
        ax.set_ylim(-0.1,1.1)
        ylab=ax.set_ylabel('Fraction')
        xlab=ax.set_xlabel('redshift')
    leg=axes[0].legend(loc=(0,1.01),ncol=1)
    leg=axes[1].legend(loc=(0,1.01),ncol=1)
    plt.savefig(fn,bbox_extra_artists=[xlab,ylab,leg], bbox_inches='tight')
    plt.close()
    print('Wrote %s' % fn)

#if __name__ == "__main__":    
grz_hist(dat)
noise_added_1(dat)
noise_added_2(dat)
delta_dec_vs_delta_ra(dat)
fraction_recovered(dat, eboss_or_desi='eboss')
e1_e2_input(dat)
e1_e2_recovered(dat)
fraction_recovered_vs_rhalf(dat)
rec_lost_contam_gr_rz_g(dat)
rec_lost_contam_grz(dat,x_ivar=0)
rec_lost_contam_by_type(dat,band='g',x_ivar=0)
rec_lost_contam_delta(dat,x_ivar=0,y_ivar=0,percentile_lines=False)
rec_lost_contam_delta_by_type(dat,band='g',
                              x_ivar=0,y_ivar=0,percentile_lines=False)
rec_lost_contam_input_elg_notelg(dat)
rec_lost_contam_fraction(dat)
    
