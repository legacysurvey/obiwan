"""Analyze the output of decals_simulations.

EXAMPLE
=======
8 500 star chunks for brick 2523p355 are here 
/project/projectdirs/desi/image_sims/2523p355
you can analyze them like this:
export DECALS_SIM_DIR=/project/projectdirs/desi/image_sims 
python legacyanalysis/decals_sim_plots.py -b 2523p355 -o STAR -out your/relative/output/path
out is optional, default is brickname/objtype

Missing object and annotated coadd plots
========================================
python legacyanalysis/decals_sim_plots.py ... --extra_plots
default is to NOT make them because chunks > 50
"""
from __future__ import division, print_function

import matplotlib
matplotlib.use('Agg') # display backend
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import matplotlib.image as mpimg

import os
import sys
import pdb
import logging
from glob import glob
import numpy as np

# import seaborn as sns
from PIL import Image, ImageDraw

def flux2mag(nmgy):
    return -2.5 * (np.log10(nmgy) - 9)

def bin_up(data_bin_by,data_for_percentile, bin_minmax=(18.,26.),nbins=20):
    '''bins "data_for_percentile" into "nbins" using "data_bin_by" to decide how indices are assigned to bins
    returns bin center,N,q25,50,75 for each bin
    '''
    bin_edges= np.linspace(bin_minmax[0],bin_minmax[1],num= nbins+1)
    vals={}
    for key in ['q50','q25','q75','n']: vals[key]=np.zeros(nbins)+np.nan
    vals['binc']= (bin_edges[1:]+bin_edges[:-1])/2.
    for i,low,hi in zip(range(nbins), bin_edges[:-1],bin_edges[1:]):
        keep= np.all((low < data_bin_by,data_bin_by <= hi),axis=0)
        if np.where(keep)[0].size > 0:
            vals['n'][i]= np.where(keep)[0].size
            vals['q25'][i]= np.percentile(data_for_percentile[keep],q=25)
            vals['q50'][i]= np.percentile(data_for_percentile[keep],q=50)
            vals['q75'][i]= np.percentile(data_for_percentile[keep],q=75)
        else:
            vals['n'][i]=0 
    return vals


def plot_flux_residual(dat,fn='flux_residual.png'):
    col = ['b', 'k', 'c', 'm', 'y', 0.8]
    fig, axes = plt.subplots(3,1, figsize=(6,8))
    plt.subplots_adjust(left=0.18,hspace=0.1)

    recovered= dat.obiwan_mask == 1
    for ax,band in zip(axes,'grz'):
        x= flux2mag(dat.get(band+'flux')[i]/\
                     dat.get('mw_transmission_'+band)[i])
        y= (dat.get('tractor_flux_'+band)[i] - dat.get(band+'flux')[i])*\
           np.sqrt(dat.get('tractor_flux_ivar_'+band)[i])
        ax.scatter(x,y,label=band,
                   s=10,edgecolor='b',c='none',lw=1.,alpha=0.5)
        ax.axhline(y=0.0,lw=1,ls='dashed',color='k')
        ax.legend(loc='upper right')
        ylab=ax.set_ylabel(r'$(\Delta Flux)/\sigma$ (Tractor - Truth)')
        xlab=ax.set_xlabel('%s mag (truth)' % band)

        binned= bin_up(x,y, bin_minmax=(0,x.max()),nbins=20)
        for perc in ['q25','q50','q75']:
            ax.plot(binned['binc'],binned[perc],c='r') 
    
    for ax in axes:
        ax.set_ylim(-10,10)
    plt.savefig(fn,bbox_extra_artists=[xlab,ylab], bbox_inches='tight')
    plt.close()
 

def basic_cut(tractor):
    '''return boolean indices for which to keep and throw out'''
    b_good= ((tractor.flux_g > 0) &
            (tractor.flux_r > 0) &
            (tractor.flux_z > 0) &
            (tractor.allmask_g == 0) &
             (tractor.allmask_r == 0) &
            (tractor.allmask_z == 0))
    b_bad= b_good == False
    return b_good, b_bad

def bright_dmag_cut(matched_simcat,matched_tractor,
                    cut=True):
    """
    Args:
        cut: True to perform good source and very bright source cuts

    Returns: 
        1) median dmag of bright sources
        2) indices of sources that are bright AND |dmag[band]| -|median_dmag[band]| > 0.005 in ANY band
        3) indices of sources that are bright AND |dmag[band]| -|median_dmag[band]| <= 0.005 in ALL bands
    """
    bright= dict(g=20.,r=19.,z=18.)
    b_good,junk= basic_cut(matched_tractor)
    # Store median dmag of bright sources
    med,b_bright={},{}
    for band in ['g','r','z']:
        # Cut to bright and good sources
        inputflux = matched_simcat.get(band+'flux')
        inputmag = 22.5-2.5*np.log10(inputflux)
        if cut:
            b_bright[band]= inputmag < bright[band]
            b_bright[band]= np.all((b_bright[band],b_good),axis=0)
        else:
            b_bright[band]= np.ones(len(matched_tractor),bool)
        #i_bright= np.where(b_bright)[0]
        #print('type(i_bright)= ',type(i_bright),"i_bright=",i_bright)
        # Compute median for each band
        inputflux = matched_simcat.get(band+'flux')[b_bright[band]]
        tractorflux = matched_tractor.get('flux_%s' % band)[b_bright[band]]
        mag_diff= -2.5*np.log10(tractorflux/inputflux)
        med[band]= np.percentile(mag_diff,q=50)
    # Boolean mask for each band
    b={}
    for band in ['g','r','z']:
        inputflux = matched_simcat.get(band+'flux')
        tractorflux = matched_tractor.get('flux_%s' % band)
        mag_diff= -2.5*np.log10(tractorflux/inputflux)
        b[band]= np.abs(mag_diff) - abs(med[band]) > 0.001
    # total boolean mask
    b_bright= np.any((b_bright['g'],b_bright['r'],b_bright['z']),axis=0)
    b_large_dmag= np.all((b_bright,b_good,b['g'],b['r'],b['z']),axis=0)
    b_small_dmag= np.all((b_bright,b_good,b['g']==False,b['r']==False,b['z']==False),axis=0)
    return med,np.where(b_large_dmag)[0],np.where(b_small_dmag)[0]
 
def plot_cutouts_by_index(simcat,index, jpeg_fn='simscoadd.jpg',
                          qafile='test.png'):
    hw = 30 # half-width [pixels]
    rad = 14
    ncols = 5
    nrows = 5
    nthumb = ncols*nrows
    dims = (ncols*hw*2,nrows*hw*2)
    mosaic = Image.new('RGB',dims)

    xpos, ypos = np.meshgrid(np.arange(0, dims[0], hw*2, dtype='int'),
                             np.arange(0, dims[1], hw*2, dtype='int'))
    im = Image.open( jpeg_fn )
    sz = im.size
    iobj = 0
    for ic in range(ncols):
        if iobj >= len(index) or iobj >= ncols*nrows: break
        for ir in range(nrows):
            if iobj >= len(index) or iobj >= ncols*nrows: break
            xx = int(simcat.x[index[iobj]])
            yy = int(sz[1]-simcat.y[index[iobj]])
            crop = (xx-hw, yy-hw, xx+hw, yy+hw)
            box = (xpos[ir, ic], ypos[ir, ic])
            thumb = im.crop(crop)
            mosaic.paste(thumb, box)
            iobj+= 1

    # Add a border and circle the missing source.
    draw = ImageDraw.Draw(mosaic)
    sz = mosaic.size
    for ic in range(ncols):
        for ir in range(nrows):
            draw.rectangle([(xpos[ir, ic], ypos[ir, ic]),
                            (xpos[ir, ic]+hw*2, ypos[ir, ic]+hw*2)])
            xx = xpos[ir, ic] + hw
            yy = ypos[ir, ic] + hw
            draw.ellipse((xx-rad, sz[1]-yy-rad, xx+rad, sz[1]-yy+rad), outline='yellow')
    mosaic.save(qafile)


def plot_annotated_coadds(simcat, jpeg_fn='simscoadd.jpg',
                          qafile='test.png'):
    rad = 7/0.262
    #imfile = os.path.join(cdir, 'qa-{}-{}-{}-{}.jpg'.format(brickname, lobjtype, suffix, chunksuffix))
    # HARDCODED fix this!!!!!
    #imfile = os.path.join(indir, 'qa-{}-{}-{}-{}.jpg'.format(brickname, lobjtype, img_name, chunksuffix))
    im = Image.open(jpeg_fn)
    sz = im.size
    draw = ImageDraw.Draw(im)
    [draw.ellipse((cat.x - rad, sz[1] - cat.y - rad,
                   cat.x + rad, sz[1] - cat.y + rad),
                  outline='yellow')
     for cat in simcat]
    im.save(qafile)


def plot_injected_mags(allsimcat, log,qafile='test.png'):
    gr_sim = -2.5*np.log10(allsimcat.gflux/allsimcat.rflux)
    rz_sim = -2.5*np.log10(allsimcat.rflux/allsimcat.zflux)
    grrange = (-0.2, 2.0)
    rzrange = (-0.4, 2.5)
    fig, ax = plt.subplots(2,1,figsize=(6,8))
    ax[0].hist(22.5-2.5*np.log10(allsimcat.rflux),bins=20,align='mid')
    ax[1].scatter(rz_sim,gr_sim,
                   s=10,edgecolor='b',c='none',lw=1.)
    for i,x_lab,y_lab in zip(range(2),['r AB','r-z'],['N','g-r']):
        xlab=ax[i].set_xlabel(x_lab)
        ylab=ax[i].set_ylabel(y_lab)
    ax[1].set_xlim(rzrange)
    ax[1].set_ylim(grrange)
    fig.subplots_adjust(wspace=0.25)
    log.info('Writing {}'.format(qafile))
    plt.savefig(qafile,bbox_extra_artists=[xlab,ylab], bbox_inches='tight')
    plt.close()
 
def plot_good_bad_ugly(bigsimcat,bigsimcat_missing, nmagbin,rminmax, log,qafile='test.png'):
    #rmaghist, magbins = np.hist a
    bigsimcat_R= flux2mag(bigsimcat.rflux)
    bigsimcat_miss_R= flux2mag(bigsimcat_missing.rflux)
    found=dict(good={},bad={},missed={})
    # bin on true r mag of matched objects, count objects in each bin
    found['rec']= bin_up(bigsimcat_R,bigsimcat_R, bin_minmax=rminmax,nbins=nmagbin) # bin_edges=magbins)
    found['los']= bin_up(bigsimcat_miss_R,bigsimcat_miss_R, bin_minmax=rminmax,nbins=nmagbin) #bin_edges=magbins)
    fig, ax = plt.subplots(1, figsize=(8,6))
    for name,color in zip(['rec','los'],['b','r']):
        ax.step(found[name]['binc'],found[name]['n'], c=color,lw=2,label=name)
    xlab=ax.set_xlabel('Input r AB')
    ylab=ax.set_ylabel('Number Recovered')
    leg=ax.legend(loc=(0.,1.01),ncol=3)
    #fig.subplots_adjust(bottom=0.15)
    log.info('Writing {}'.format(qafile))
    plt.savefig(qafile, bbox_extra_artists=[leg,xlab,ylab], bbox_inches='tight')
    plt.close()

def plot_tractor_minus_answer(bigsimcat,bigtractor, rminmax, log,qafile='test.png'):
    b_good,junk= basic_cut(bigtractor)    
    fig, ax = plt.subplots(3, sharex=True, figsize=(6,8))

    col = ['b', 'k', 'c', 'm', 'y', 0.8]
    rmag = bigsimcat.rflux
    for thisax, thiscolor, band in zip(ax, col, ('g','r','z')):
        inputflux = bigsimcat.get(band+'flux')
        tractorflux = bigtractor.get('flux_%s' % band)
        tractorivar = bigtractor.get('flux_ivar_%s' % band)
        #import pickle
        #fout=open('test.pickle','w')
        #pickle.dump((tractorflux,inputflux,b_good),fout)
        #fout.close()
        #print('exiting early')
        #sys.exit()
        inputmag = 22.5-2.5*np.log10(inputflux[b_good])
        mag_diff= -2.5*np.log10(tractorflux[b_good]/inputflux[b_good])
        thisax.scatter(inputmag, mag_diff,
                       s=10,edgecolor=thiscolor,c='none',lw=1.)
        thisax.set_ylim(-0.1,0.1)
        thisax.set_xlim(inputmag.min()-0.1, inputmag.max()+0.1)
        thisax.axhline(y=0.0,lw=2,ls='solid',color='gray')
        #arr= np.ma.masked_array(mag_diff, mask= np.isfinite(mag_diff) == False)
        # Compute median
        inputflux = bigsimcat.get(band+'flux')
        tractorflux = bigtractor.get('flux_%s' % band)
        mag_diff= -2.5*np.log10(tractorflux/inputflux)
        med= np.percentile(mag_diff,q=50)
        
        thisax.axhline(y=med,lw=2,ls='dashed',color='red',label='Median=%.3f' % med)
        thisax.legend(loc='upper left',fontsize='x-small')
        #thisax.text(0.05,0.05, band.lower(), horizontalalignment='left',
                    #verticalalignment='bottom',transform=thisax.transAxes,
                    #fontsize=16)
    ax[0].set_ylabel('$\Delta$g')
    ax[1].set_ylabel('$\Delta$r (Tractor - Input)')
    ylab=ax[2].set_ylabel('$\Delta$z')
    xlab=ax[2].set_xlabel('Input magnitude (AB mag)')
    fig.subplots_adjust(left=0.18,hspace=0.1)
    log.info('Writing {}'.format(qafile))
    plt.savefig(qafile,bbox_extra_artists=[xlab,ylab], bbox_inches='tight')
    plt.close()


def plot_color_tractor_minus_answer(bigsimcat, bigtractor,rminmax, brickname,lobjtype, log,qafile='test.png'):
    gr_tra = -2.5*np.log10(bigtractor.flux_g/bigtractor.flux_r)
    rz_tra = -2.5*np.log10(bigtractor.flux_r/bigtractor.flux_z)
    gr_sim = -2.5*np.log10(bigsimcat.gflux/bigsimcat.rflux)
    rz_sim = -2.5*np.log10(bigsimcat.rflux/bigsimcat.zflux)
    rmag = flux2mag(bigsimcat.rflux)

    col = ['b', 'k', 'c', 'm', 'y', 0.8]
    fig, ax = plt.subplots(2,sharex=True,figsize=(6,8))
    
    ax[0].scatter(rmag, gr_tra-gr_sim, color=col[0], s=10)
    ax[1].scatter(rmag, rz_tra-rz_sim, color=col[1], s=10)
    
    [thisax.set_ylim(-0.7,0.7) for thisax in ax]
    [thisax.set_xlim(rminmax + [-0.1, 0.0]) for thisax in ax]
    [thisax.axhline(y=0.0, lw=2, ls='solid', color='gray') for thisax in ax]
    
    ax[0].set_ylabel('$\Delta$(g - r) (Tractor minus Input)')
    ax[1].set_ylabel('$\Delta$(r - z) (Tractor minus Input)')
    ax[1].set_xlabel('Input r magnitude (AB mag)')
    fig.subplots_adjust(left=0.18,hspace=0.1)

    log.info('Writing {}'.format(qafile))
    plt.savefig(qafile)
    plt.close()

def plot_fraction_recovered(allsimcat,bigsimcat, nmagbin,rminmax, brickname, lobjtype, log,qafile='test.png'):
    allsimcat_R= flux2mag(allsimcat.rflux)
    bigsimcat_R= flux2mag(bigsimcat.rflux)
    rmaghist, magbins = np.histogram(allsimcat_R, bins=nmagbin, range=rminmax)
    cmagbins = (magbins[:-1] + magbins[1:]) / 2.0
    ymatch, binsmatch = np.histogram(bigsimcat_R, bins=nmagbin, range=rminmax)
    fig, ax = plt.subplots(1, figsize=(8,6))
    ax.step(cmagbins, 1.0*ymatch/rmaghist, c='k',lw=3,label='All objects')
    #ax.step(cmagbins, 1.0*ymatchgood/rmaghist, lw=3, ls='dashed', label='|$\Delta$m|<0.3')
    ax.axhline(y=1.0,lw=2,ls='dashed',color='k')
    ax.set_xlabel('Input r magnitude (AB mag)')
    ax.set_ylabel('Fraction Recovered'.format(lobjtype))
    ax.set_ylim([0.0, 1.1])
    ax.legend('lower left')
    fig.subplots_adjust(bottom=0.15)
    log.info('Writing {}'.format(qafile))
    plt.savefig(qafile)
    plt.close()

def plot_sn_recovered(allsimcat,bigsimcat,bigtractor, brickname, lobjtype, log,qafile='test.png'):
    allsimcat_R= flux2mag(allsimcat.rflux)
    # min,max mag of all bands
    grrange = (-0.2, 2.0)
    rzrange = (-0.4, 2.5)
    rmin,rmax= allsimcat_R.min(), allsimcat_R.max()
    mag_min= np.min((rmin,rmin+grrange[0],rmin-rzrange[1]))
    mag_max= np.max((rmax,rmax+grrange[1],rmax-rzrange[0]))
    s2n=dict(g={},r={},z={})
    for band in ['g','r','z']:
        mag= 22.5-2.5*np.log10(bigsimcat.get(band+'flux'))
        # HARDCODED mag range
        s2n[band]= bin_up(mag, bigtractor.get('flux_%s' % band)*np.sqrt(bigtractor.get('flux_ivar_%s' % band)), \
                          bin_minmax=(18,26),nbins=20) 
    fig, ax = plt.subplots(1, figsize=(8,6))
    xlab=ax.set_xlabel('Input magnitude (AB)')
    ylab=ax.set_ylabel(r'Median S/N = $F/\sigma$',fontweight='bold',fontsize='large')
    title= ax.set_title('S/N of Recovered Objects')
    for band,color in zip(['g','r','z'],['g','r','b']):
        ax.plot(s2n[band]['binc'], s2n[band]['q50'],c=color,ls='-',lw=2,label=band)
        #ax.fill_between(s2n[band]['cbin'],s2n[band]['q25'],s2n[band]['q75'],color=color,alpha=0.25)
    ax.axhline(y=5.,lw=2,ls='dashed',color='k',label='S/N = 5')
    ax.set_yscale('log')
    leg=ax.legend(loc=3)
    fig.subplots_adjust(bottom=0.15)
    log.info('Writing {}'.format(qafile))
    plt.savefig(qafile,bbox_extra_artists=[leg,xlab,ylab,title], bbox_inches='tight',dpi=150)
    plt.close()

def plot_recovered_types(bigsimcat,bigtractor, nmagbin,rminmax, objtype,log,qafile='test.png'):
    bigsimcat_R= flux2mag(bigsimcat.rflux)
    fig = plt.figure(figsize=(8, 6))
    ax = fig.gca()
    rmaghist, magbins = np.histogram(bigsimcat_R, bins=nmagbin, range=rminmax)
    cmagbins = (magbins[:-1] + magbins[1:]) / 2.0
    tractortype = np.char.strip(bigtractor.get('type'))
    for otype in ['PSF', 'SIMP', 'EXP', 'DEV', 'COMP']:
        these = np.where(tractortype == otype)[0]
        if len(these)>0:
            yobj, binsobj = np.histogram(bigsimcat_R[these], bins=nmagbin, range=rminmax)
            #plt.step(cmagbins,1.0*yobj,lw=3,alpha=0.5,label=otype)
            plt.step(cmagbins,1.0*yobj/rmaghist,lw=3,alpha=0.5,label=otype)
    plt.axhline(y=1.0,lw=2,ls='dashed',color='gray')
    plt.xlabel('Input r magnitude (AB mag)')
    #plt.ylabel('Number of Objects')
    plt.ylabel('Fraction of {}s classified'.format(objtype))
    plt.ylim([0.0,1.1])
    plt.legend(loc='center left', bbox_to_anchor=(0.08,0.5))
    fig.subplots_adjust(bottom=0.15)

    log.info('Writing {}'.format(qafile))
    plt.savefig(qafile)
    plt.close()


 
def create_confusion_matrix(answer_type,predict_type, types=['PSF','SIMP','EXP','DEV','COMP','REX'],slim=True):
    '''compares classifications of matched objects, returns 2D array which is conf matrix and xylabels
    return 5x5 confusion matrix and colum/row names
    answer_type,predict_type -- arrays of same length with reference and prediction types'''
    for typ in set(answer_type): assert(typ in types)
    for typ in set(predict_type): assert(typ in types)
    # if a type was not in answer (training) list then don't put in cm
    if slim: ans_types= set(answer_type)
    # put in cm regardless
    else: ans_types= set(types)
    cm=np.zeros((len(ans_types),len(types)))-1
    for i_ans,ans_type in enumerate(ans_types):
        ind= np.where(answer_type == ans_type)[0]
        for i_pred,pred_type in enumerate(types):
            n_pred= np.where(predict_type[ind] == pred_type)[0].size
            if ind.size > 0: cm[i_ans,i_pred]= float(n_pred)/ind.size # ind.size is constant for loop over pred_types
            else: cm[i_ans,i_pred]= np.nan
    if slim: return cm,ans_types,types #size ans_types != types
    else: return cm,types

def plot_confusion_matrix(cm,answer_names,all_names, log, qafile='test.png'):
    plt.imshow(cm, interpolation='nearest', cmap=plt.cm.Blues, vmin=0,vmax=1)
    cbar=plt.colorbar()
    plt.xticks(range(len(all_names)), all_names)
    plt.yticks(range(len(answer_names)), answer_names)
    ylab=plt.ylabel('True')
    xlab=plt.xlabel('Predicted (tractor)')
    for row in range(len(answer_names)):
        for col in range(len(all_names)):
            if np.isnan(cm[row,col]): 
                plt.text(col,row,'n/a',va='center',ha='center')
            elif cm[row,col] > 0.5: 
                plt.text(col,row,'%.2f' % cm[row,col],va='center',ha='center',color='yellow')
            else: 
                plt.text(col,row,'%.2f' % cm[row,col],va='center',ha='center',color='black')
    log.info('Writing {}'.format(qafile))
    plt.savefig(qafile, bbox_extra_artists=[xlab,ylab], bbox_inches='tight',dpi=150)
    plt.close()

def plot_cm_stack(cm_stack,stack_names,all_names, log, qafile='test.png'):
    '''cm_stack -- list of single row confusion matrices
    stack_names -- list of same len as cm_stack, names for each row of cm_stack'''
    # combine list into single cm
    cm=np.zeros((len(cm_stack),len(all_names)))+np.nan
    for i in range(cm.shape[0]):
        if len(cm_stack[i]) > 0:
            cm[i,:]= cm_stack[i]
    # make usual cm, but labels repositioned
    plt.imshow(cm, interpolation='nearest', cmap=plt.cm.Blues)
    cbar=plt.colorbar()
    plt.xticks(range(len(all_names)), all_names)
    plt.yticks(range(len(stack_names)), stack_names)
    ylab=plt.ylabel('True = PSF')
    xlab=plt.xlabel('Predicted (tractor)')
    for row in range(len(stack_names)):
        for col in range(len(all_names)):
            if np.isnan(cm[row,col]): 
                plt.text(col,row,'n/a',va='center',ha='center')
            elif cm[row,col] > 0.5: 
                plt.text(col,row,'%.2f' % cm[row,col],va='center',ha='center',color='yellow')
            else: 
                plt.text(col,row,'%.2f' % cm[row,col],va='center',ha='center',color='black')
            #if np.isnan(cm[row,col]): 
            #    plt.text(col,row,'n/a',va='center',ha='center')
            #else: plt.text(col,row,'%.2f' % cm[row,col],va='center',ha='center')
    log.info('Writing {}'.format(qafile))
    plt.savefig(qafile, bbox_extra_artists=[xlab,ylab], bbox_inches='tight',dpi=150)
    plt.close()

def make_stacked_cm(bigsimcat,bigtractor, log,qafile='test.png'):
    bigsimcat_R= flux2mag(bigsimcat.rflux)
    types= ['PSF', 'REX','SIMP', 'EXP', 'DEV', 'COMP']
    cm_stack,stack_names=[],[]
    rbins= np.array([18.,20.,22.,23.,24.])
    for rmin,rmax in zip(rbins[:-1],rbins[1:]):
        # master cut
        br_cut= ((bigsimcat_R > rmin) &
                 (bigsimcat_R <= rmax))
        stack_names+= ["%d < r <= %d" % (int(rmin),int(rmax))]
        cm,ans_names,all_names= create_confusion_matrix(np.array(['PSF']*len(bigtractor))[br_cut],
                                                        np.char.strip(bigtractor.get('type'))[br_cut], \
                                                        types=types)
        cm_stack+= [cm]
    plot_cm_stack(cm_stack, stack_names,all_names, log, qafile=qafile)

class TestData(object):
    def fetch(self,download_dir,outdir):
        name='qa.tar.gz'
        remote_fn= os.path.join(download_dir,name)
        fetch_targz(remote_fn,outdir)

def simcatfn_to_rs_dir(simcat_fn):
    """returns abs path to rs* or skip_rs* directory given full path to simcat_fn"""
    return os.path.dirname( os.path.dirname( simcat_fn))
    
class Checks(object):
    """various sanity checks on obiwan outputs    
    """
    def __init__(self,warn=False):
        """
        Args:
            warn: raise warning instead of error
        """
        self.warn= warn
    
    def test_nobj(self,metacat_fn,simcat_fn,skipid_fn):
        meta= fits_table(metacat_fn)
        simcat= fits_table(simcat_fn)
        if os.path.exists(skipid_fn):
            skipid= fits_table(skipid_fn)
        else:
            skipid= []
        print(metacat_fn,simcat_fn,skipid_fn)
        print(meta.nobj,len(simcat),len(skipid))
        #if self.warn:
        #    if meta.nobj != len(simcat) + len(skipid):
        #        print('WARNING: %d %d %d' %
        #              (meta.nobj,len(simcat),len(skipid)))
        #else:
        # simcat+skipid < 300 if ran out of randoms, 300 otherwise
        assert(len(simcat) + len(skipid) <= meta.nobj)


class SourceMatcher(object):
    """Does all the matching between injected, recovered, and pre existing sources
    
    There are 6 types of sources. 
        (1) added and recovered 
        (2) added but lost
        (3) pre existng from previous DR and recovered even when simulated sources injected, 
        (4) '' but lost for whatever reason,  
        (5) added and pre existing near by and found at least one of them
        (6) '' but lost both
    """
    def __init__(self,verbose=True):
        """
        Attrs:
            verbose: True to print stats as go along
            i_catalogue_name: indices for that catalogue
            catalouge_name: that catalogue with those indices applied
            size: dict number sources in each tyep fo catalogue
        """
        self.verbose= verbose
        
        self.i_simcat= {}
        self.simcat= {}
        
        self.i_simtractor= {}
        self.simtractor= {}

        self.i_realtractor= {}
        self.realtractor= {}

        self.size= {}

    def added(self,simcat,simtractor):
        """Source types 1, 2
        
        Args:
            simcat: all injected sources, before geometry cuts
            simtractor: all sources in resulting tractor catalogues
        """
        self.size['simcat']= len(simcat)
        self.size['simtractor']= len(simtractor)

        I,J,d = match_radec(simcat.ra, simcat.dec,
                            simtractor.ra, simtractor.dec,
                            1./3600, nearest=True)
        self.i_simtractor['rec_simcat']= J #indices of simtractor, where recovered simcat sources
        self.i_simcat['recby_simtractor']= I #indices of simcat, where recovered by simtractor
        i_simcat_all= np.arange(self.size['simcat'])
        self.i_simcat['losby_simtractor']= list( set(i_simcat_all).difference(set(I)) )

        # Apply
        for key in self.i_simcat.keys():
            self.simcat[key]= simcat.copy()[ self.i_simcat[key] ]
        for key in self.i_simtractor.keys():
            self.simtractor[key]= simtractor.copy()[ self.i_simtractor[key] ]
        if self.verbose:
            self.print_added()

    def print_added(self):
        assert( hasattr(self,'i_simtractor') )
        print('Added and recovered: %d/%d, Added but lost: %d/%d' % 
               (len(self.i_simcat['recby_simtractor']), self.size['simcat'],
                len(self.i_simcat['losby_simtractor']), self.size['simcat']))

    def already_exist(self,realtractor,simtractor):
        """Source types 3, 4
        
        WARNING: assumes the obiwan run used identical in put ccds and code to 
            the data release that realtractor comes from

        Args:
            realtractor: tractor catalogue from a data release
            simtractor: all sources in resulting tractor catalogues
        """
        self.size['realtractor']= len(realtractor)
        self.size['simtractor']= len(simtractor)

        I,J,d = match_radec(realtractor.ra, realtractor.dec,
                            simtractor.ra, simtractor.dec,
                            1./3600, nearest=True)

        self.i_simtractor['rec_butreal']= J 
        self.i_realtractor['recby_simtractor']= I 
        i_realtractor_all= np.arange(self.size['realtractor'])
        self.i_realtractor['losby_simtractor']= list( set(i_realtractor_all).difference(set(I)) )

        # Apply
        for key in self.i_realtractor.keys():
            self.realtractor[key]= realtractor.copy()[ self.i_realtractor[key] ]
        for key in self.i_simtractor.keys():
            self.simtractor[key]= simtractor.copy()[ self.i_simtractor[key] ]
        if self.verbose:
            self.print_real()

    def print_real(self):
        assert( hasattr(self,'i_realtractor') )
        print('Real and recovered: %d/%d, Real but lost: %d/%d' % 
                (len(self.i_realtractor['recby_simtractor']), self.size['realtractor'],
                 len(self.i_realtractor['losby_simtractor']), self.size['realtractor']))
    
    # Added and Existing source there
    # I,J,d = match_radec(pre_exist.ra,pre_exist.dec,
    #                     simcat.ra,simcat.dec,
    #                     1./3600, nearest=True)
    # i_add_exist_all= np.arange(len(I))
    # i_add_exist_fnd= list( set(J).intersection(set(i_add_all)) )
    # i_exist_lost= list( set(i_exist_all).difference(set(i_exist_fnd)) )
    # exist_fnd= simtractor.copy().cut(i_exist_fnd)
    # exist_lost= simtractor.copy().cut(i_exist_lost)
    # print('Exists and found: %d/%d, Exists but lost: %d/%d' % 
    #       (len(exist_fnd),len(pre_exist),len(exist_lost),len(pre_exist)))
  
def imshow_one_chunk(simcat_fn, simtractor_fn, tractor_fn):
    simcat= fits_table(simcat_fn)
    simtractor= fits_table(simtractor_fn)
    data= SourceMatcher()
    data.added(simcat,simtractor)
    
    chunk_name= os.path.basename( simcatfn_to_rs_dir(simcat_fn))
    # Circle added
    for img_name in ('simscoadd','image', 'resid'):
        jpeg_fn= os.path.join(simcatfn_to_rs_dir(simcat_fn),
                              'coadd/legacysurvey-%s-%s.jpg' %
                              (brickname,img_name))
        qafile = os.path.join(output_dir, 
                              'qa_{}_{}_{}_{}_annot.png'.format(brickname, objtype,img_name,chunk_name))
        plot_annotated_coadds(simcat, jpeg_fn=jpeg_fn, qafile=qafile)    
        log.info('Wrote {}'.format(qafile))
    
    # Stamps added and recovered
    for img_name in ['image','resid','simscoadd']:
        jpeg_fn= os.path.join(simcatfn_to_rs_dir(simcat_fn),
                              'coadd/legacysurvey-%s-%s.jpg' %
                              (brickname,img_name))
        qafile = os.path.join(output_dir, 
                              'qa_{}_{}_{}_bright_large_dmag_{}.png'.format(brickname, objtype, img_name,chunk_name))
        # Large dmag cutouts
        junk,i_large_dmag,i_small_dmag= (
            bright_dmag_cut(data.simcat['recby_simtractor'],
                            data.simtractor['rec_simcat'],
                            cut=False) )        
        plot_cutouts_by_index(simcat,i_large_dmag,
                              jpeg_fn=jpeg_fn, qafile=qafile)
        log.info('Wrote {}'.format(qafile))
        # Small dmag cutouts
        qafile = os.path.join(output_dir, 
                              'qa_{}_{}_{}_bright_small_dmag_{}.png'.format(brickname, objtype, img_name,chunk_name))
        plot_cutouts_by_index(simcat,i_small_dmag,
                              jpeg_fn=jpeg_fn, qafile=qafile)
        log.info('Wrote {}'.format(qafile))
        
    # Stamps added but lost
    if len( data.i_simcat['losby_simtractor']) > 0:
        simcat_R= flux2mag(simcat.rflux)
        for img_name in ['image']: #,'simscoadd']:
            jpeg_fn= os.path.join(simcatfn_to_rs_dir(simcat_fn),
                                  'coadd/legacysurvey-%s-%s.jpg' %
                                  (brickname,img_name))
            qafile = os.path.join(output_dir, 
                                  'qa_{}_{}_{}_missing_{}.png'.format(brickname, objtype, img_name,chunk_name))
            miss = data.i_simcat['losby_simtractor']
            plot_cutouts_by_index(simcat,miss,
                                  jpeg_fn=jpeg_fn, qafile=qafile)
            log.info('Wrote {}'.format(qafile))
        

if __name__ == "__main__":    
    from argparse import ArgumentParser
    parser = ArgumentParser(description='DECaLS simulations.')
    parser.add_argument('--data_dir', default='~/mydata', help='where obiwan outputs and DR5 outputs have been untarred to',required=False)
    parser.add_argument('--brick',help='process this brick',required=True)
    parser.add_argument('--objtype', type=str, choices=['star','qso','elg','lrg'], help='object type',required=True) 
    parser.add_argument('--output_dir', default='./',help='relative path to output directory',required=False) 
    parser.add_argument('--skip_coadds', action='store_true',default=False, 
                        help='speed up by not drawing added sources on coadd jpegs')
    parser.add_argument('--verbose', action='store_true', 
                        help='toggle on verbose output')
    args = parser.parse_args()

    from obiwan.fetch import fetch_targz
    DOWNLOAD_ROOT = "http://portal.nersc.gov/project/desi/users/kburleigh/obiwan/"
    Tdata= TestData()
    Tdata.fetch(DOWNLOAD_ROOT,args.data_dir)
    
    # Set the debugging level
    if args.verbose:
        lvl = logging.DEBUG
    else:
        lvl = logging.INFO
    logging.basicConfig(format='%(message)s', level=lvl, stream=sys.stdout)
    log = logging.getLogger('__name__')

    brickname = args.brick
    objtype = args.objtype
    log.info('Analyzing objtype {} on brick {}'.format(objtype, brickname))

    kwargs= dict(do_skipids= 'no',
                 decals_sim_dir= os.path.join(args.data_dir,
                                              'qa','elg_9deg2_ra175'),
                 objtype= args.objtype,
                 brickname= args.brick,
                 rowst=7)
    input_dir= get_savedir(**kwargs)
    input_dir= os.path.dirname(input_dir)

    output_dir= os.path.join(input_dir, 'qaplots')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    #metacat_fns= glob( os.path.join(input_dir,
    #                               '*/obiwan/metacat-%s-%s.fits' %
    #                               (args.objtype, args.brick)))
    simcat_fns= glob( os.path.join(input_dir,
                                   '*/obiwan/simcat-%s-%s.fits' %
                                   (args.objtype, args.brick)))
    #simtractor_fns= glob( os.path.join(input_dir,
    #                                   '*/tractor/tractor-%s.fits' %
    #                                   (args.brick,)))
    realtractor_fn= os.path.join(args.data_dir,
                             'qa/DR5_out/tractor/%s/tractor-%s.fits' %
                             (args.brick[:3],args.brick) )
    assert( (len(simcat_fns) > 0) &
            (len(realtractor_fn) > 0) )
    metacat_fns= [fn.replace('simcat','metacat')
                  for fn in simcat_fns]
    skipid_fns= [fn.replace('simcat','skippedids')
                 for fn in simcat_fns]
    simtractor_fns= [os.path.dirname(fn).replace(
                         'obiwan','tractor/tractor-%s.fits' %
                         args.brick)
                     for fn in simcat_fns]

    # Sanity checks
    sanity= Checks(warn=False)
    for metacat_fn, simcat_fn, skipid_fn in zip(
            metacat_fns, simcat_fns, skipid_fns):
      sanity.test_nobj(metacat_fn,simcat_fn, skipid_fn)   

    # Plots for each rs*, skip_rs*
    if not args.skip_coadds:
        for simcat_fn, simtractor_fn in zip(
                simcat_fns,simtractor_fns):
            imshow_one_chunk(simcat_fn, simtractor_fn, 
                             realtractor_fn)

    # Merge
    simcat= stack_tables(simcat_fns,textfile=False)
    simtractor= stack_tables(simtractor_fns,textfile=False)
    realtractor= fits_table(realtractor_fn)
    print('Injected %d sources, Real %d sources, Recovered sims + real %d' % 
          (len(simcat),len(realtractor),len(simtractor)))

    data= SourceMatcher()
    data.added(simcat,simtractor)
    data.already_exist(realtractor,simtractor)

    ###############
    # Tally for all rs*, skip_rs*
    temp=EmptyClass()
    temp.realtractor= fits_table(realtractor_fn)
    print('-'*40)
    print('%s %10s %s/%s %s/%s' %
          (brickname,'chunk','rec','inj','rec','real'))
    for simcat_fn, simtractor_fn in zip(
            simcat_fns,simtractor_fns):
        temp.simcat= fits_table(simcat_fn)
        temp.simtractor= fits_table(simtractor_fn)

        mat= SourceMatcher(verbose=False)
        mat.added(temp.simcat,temp.simtractor)
        mat.already_exist(temp.realtractor,temp.simtractor)

        rs= os.path.basename(simcatfn_to_rs_dir(simcat_fn))
        print('\t %10s %d/%d %d/%d' %
              (rs,
               len(mat.i_simcat['recby_simtractor']), mat.size['simcat'],
               len(mat.i_realtractor['recby_simtractor']), mat.size['realtractor']))
    print('\t %10s %d/%d %d/%d' %
          ('total'.upper(),
           len(data.i_simcat['recby_simtractor']), data.size['simcat'],
           len(data.i_realtractor['recby_simtractor']), data.size['realtractor']))
    print('-'*40)
    #############
    raise ValueError
    
    
    #good = np.where((np.abs(tractor['decam_flux'][m1,2]/simcat['rflux'][m2]-1)<0.3)*1)

    # Plotting preferences
    #sns.set(style='white',font_scale=1.6,palette='dark')#,font='fantasy')
    #col = sns.color_palette('dark')
    col = ['b', 'k', 'c', 'm', 'y', 0.8]
    
    # We need this for our histograms below
    magbinsz = 0.2
    # See priors.py
    mag_limits= dict(elg_rlimit=23.4+1,
                     lrg_zlimit=20.46+1)
    rminmax = np.array([15., mag_limits['%s_rlimit' % args.objtype]])  
    nmagbin = long((rminmax[1]-rminmax[0])/magbinsz)



    # now operate on concatenated catalogues from multiple chunks
    # Grab flags
    b_good,b_bad= basic_cut(simtractor)

    # mags and colors of ALL injected sources
    plot_injected_mags(simcat, log, 
                       qafile= os.path.join(output_dir, 
                                            'qa-{}-{}-injected-mags.png'.format(
                                             brickname, objtype)))
   
    # number of recovered sources that are bad, good and number of lost, binned by r mag
    plot_good_bad_ugly(data.simcat['recby_simtractor'],
                       data.simcat['losby_simtractor'], 
                       nmagbin,rminmax, log, 
                       qafile= os.path.join(output_dir, 
                                            'qa-{}-{}-N-good-bad-missed.png'.format(
                                             brickname, objtype)))

    # Flux residuals vs input magnitude
    plot_tractor_minus_answer(data.simcat['recby_simtractor'], 
                              data.simtractor['rec_simcat'], 
                              rminmax, log, 
                              qafile= os.path.join(output_dir, 
                                                   'qa-{}-{}-good-flux.png'.format(
                                                    brickname, objtype)))
 
    # chi plots: Flux residual / estimated Flux error
    plot_chi(data.simcat['recby_simtractor'], 
             data.simtractor['rec_simcat'], 
             rminmax, log,
             qafile= os.path.join(output_dir,
                                  'qa-{}-{}-chi-good.png'.format(brickname, objtype)))
   
    # Color residuals
    plot_color_tractor_minus_answer(data.simcat['recby_simtractor'], 
                                    data.simtractor['rec_simcat'],
                                    rminmax, brickname,objtype, log, qafile =\
             os.path.join(output_dir, 'qa-{}-{}-color.png'.format(brickname, objtype)))

    # Fraction of recovered sources
    plot_fraction_recovered(simcat, 
                            data.simcat['recby_simtractor'], 
                            nmagbin,rminmax, brickname, objtype, log, qafile =\
            os.path.join(output_dir, 'qa-{}-{}-frac.png'.format(brickname, objtype)))

    # S/N of recovered sources (S/N band vs. AB mag band)
    plot_sn_recovered(simcat, 
                      data.simcat['recby_simtractor'],
                      data.simtractor['rec_simcat'],
                      brickname, objtype, log, qafile =\
             os.path.join(output_dir, 'qa-{}-{}-SN.png'.format(brickname, objtype)))

    # Distribution of object types for matching sources.
    plot_recovered_types(data.simcat['recby_simtractor'],
                      data.simtractor['rec_simcat'],
                      nmagbin,rminmax, objtype,log, qafile =\
             os.path.join(output_dir, 'qa-{}-{}-type.png'.format(brickname, objtype)))

    # Confusion matrix for distribution of object types
    # Basic cm, use slim=False
    types= ['PSF','REX', 'SIMP', 'EXP', 'DEV', 'COMP']
    cm,all_names= create_confusion_matrix(np.array(['PSF']*len(data.simtractor['rec_simcat'])),
                                          np.char.strip(data.simtractor['rec_simcat'].get('type')), 
                                          types=types,slim=False)
    qafile = os.path.join(output_dir, 'qa-{}-{}-confusion.png'.format(brickname, objtype))
    plot_confusion_matrix(cm,all_names,all_names, log,qafile)
    
    # Now a stacked confusion matrix
    # Compute a row for each r mag range and stack rows
    make_stacked_cm(data.simcat['recby_simtractor'],
                    data.simtractor['rec_simcat'],
                    log,qafile =\
             os.path.join(output_dir, 'qa-{}-{}-confusion-stack.png'.format(brickname, objtype)))
   
    '''
    # Morphology plots
    if objtype=='ELG':
        fig = plt.figure(figsize=(8,4))
        plt.subplot(1,3,1)
        plt.plot(rmag,deltam,'s',markersize=3)
        plt.axhline(y=0.0,lw=2,ls='solid',color='gray')
        plt.xlim(rminmax)
        plt.xlabel('r (AB mag)')

        plt.subplot(1,3,2)
        plt.plot(bigsimcat['R50_1'],deltam,'s',markersize=3)
        plt.axhline(y=0.0,lw=2,ls='solid',color='gray')
        plt.xlabel('$r_{50}$ (arcsec)')

        plt.subplot(1,3,3)
        plt.plot(bigsimcat['BA_1'],deltam,'s',markersize=3)
        plt.axhline(y=0.0,lw=2,ls='solid',color='gray')
        plt.xlabel('b/a')
        plt.xlim([0.2,1.0])
        fig.subplots_adjust(bottom=0.18)
        qafile = os.path.join(output_dir,'qa-'+brickname+'-'+lobjtype+'-morph.png')
        log.info('Writing {}'.format(qafile))
        plt.savefig(qafile)
    '''
