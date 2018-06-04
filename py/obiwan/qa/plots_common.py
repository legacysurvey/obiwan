import matplotlib
matplotlib.use('Agg') # display backend
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import matplotlib.image as mpimg
import numpy as np
import pandas as pd

def flux2mag(nmgy):
    """return ABmag"""
    return -2.5 * (np.log10(nmgy) - 9)

def mag2flux(ABmag):
    """return nmgy"""
    return 10**(-ABmag/2.5 + 9)

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

def myhist2D(ax,x,y,xlim=(),ylim=(),nbins=()):
    #http://www.astroml.org/book_figures/chapter1/fig_S82_hess.html
    H, xbins, ybins = np.histogram2d(x,y,
                                     bins=(np.linspace(xlim[0],xlim[1],nbins[0]),
                                           np.linspace(ylim[0],ylim[1],nbins[1])))
    # Create a black and white color map where bad data (NaNs) are white
    cmap = plt.cm.binary
    cmap.set_bad('w', 1.)
    
    H[H == 0] = 1  # prevent warnings in log10
    ax.imshow(np.log10(H).T, origin='lower',
              #extent=[xbins[0], xbins[-1], ybins[0], ybins[-1]],
              extent=[xlim[0], xlim[1], ylim[0], ylim[1]],
              cmap=cmap, interpolation='nearest',
              aspect='auto')
    #ax.set_xlim(xlim)
    #ax.set_ylim(ylim)

def myhist(ax,data,bins=20,color='b',normed=False,ls='-',label=None):
    if label:
        _=ax.hist(data,bins=bins,facecolor=color,normed=normed,
                  histtype='stepfilled',edgecolor='none', alpha=0.75,label=label)
    else:
        _=ax.hist(data,bins=bins,facecolor=color,normed=normed,
                  histtype='stepfilled',edgecolor='none', alpha=0.75)

def myhist_step(ax,data,bins=20,color='b',normed=False,lw=2,ls='solid',label=None,
                return_vals=False,domain=None):
    kw= dict(bins=bins,color=color,normed=normed,
             histtype='step',range=domain,lw=lw,ls=ls)
    if label:
        kw.update(label=label)
    h,bins,_=ax.hist(data,**kw)
    if return_vals:
        return h,bins

def myscatter(ax,x,y, color='b',m='o',s=10.,alpha=0.75,label=None):
    if label is None or label == 'None':
        ax.scatter(x,y, c=color,edgecolors='none',marker=m,s=s,rasterized=True,alpha=alpha)
    else:
        ax.scatter(x,y, c=color,edgecolors='none',marker=m,s=s,rasterized=True,alpha=alpha,label=label)

def myscatter_open(ax,x,y, color='b',m='o',s=10.,alpha=0.75,label=None):
    if label is None or label == 'None':
        ax.scatter(x,y, facecolors='none',edgecolors=color,marker=m,s=s,rasterized=True,alpha=alpha)
    else:
        ax.scatter(x,y, facecolors='none',edgecolors=color,marker=m,s=s,rasterized=True,alpha=alpha,label=label)



def myannot(ax,xarr,yarr,sarr, ha='left',va='bottom',fontsize=20):
    '''x,y,s are iterable'''
    for x,y,s in zip(xarr,yarr,sarr):
        ax.annotate(s,xy=(x,y),xytext=(x,y),
                    horizontalalignment=ha,verticalalignment=va,fontsize=fontsize)

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


def myerrorbar(ax,x,y, yerr=None,xerr=None,color='b',ls='none',m='o',s=10.,mew=1,alpha=0.75,label=None):
    if label is None or label == 'None':
        ax.errorbar(x,y, xerr=xerr,yerr=yerr,ls=ls,alpha=alpha,
                    marker=m,ms=s,mfc=color,mec=color,ecolor=color,mew=mew)
    else:
        ax.errorbar(x,y, xerr=xerr,yerr=yerr,ls=ls,alpha=alpha,label=label,
                    marker=m,ms=s,mfc=color,mec=color,ecolor=color,mew=mew)

def create_confusion_matrix(answer_type,predict_type, 
                            poss_ans_types=['PSF','SIMP','EXP','DEV','COMP','REX'],
                            poss_pred_types=['PSF','SIMP','EXP','DEV','COMP','REX']):
    '''compares classifications of matched objects, returns 2D array which is conf matrix and xylabels
    return 5x5 confusion matrix and colum/row names
    answer_type,predict_type -- arrays of same length with reference and prediction types'''
    for typ in set(answer_type): assert(typ in poss_ans_types)
    for typ in set(predict_type): 
        print('typ=',typ)
        assert(typ in poss_pred_types)
    cm=np.zeros((len(poss_ans_types),len(poss_pred_types)))-1
    for i_ans,ans_type in enumerate(poss_ans_types):
        ind= np.where(answer_type == ans_type)[0]
        for i_pred,pred_type in enumerate(poss_pred_types):
            n_pred= np.where(predict_type[ind] == pred_type)[0].size
            if ind.size > 0: cm[i_ans,i_pred]= float(n_pred)/ind.size # ind.size is constant for loop over pred_types
            else: cm[i_ans,i_pred]= 0 #np.nan
    return cm

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

def to_csv(d,fn='test.csv'):
    df= pd.DataFrame(d)
    #df= df.round({})
    df.to_csv(fn,index=False)
    print('Wrote %s' % fn)


