'''
All plotting code for theValidator
'''
if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg') #display backend
import matplotlib.pyplot as plt
import os
import sys
import numpy as np
import scipy.stats

from tractor.brightness import NanoMaggies
import statsmodels.api as sm
from scipy.optimize import curve_fit

from catalogues import Cuts4MatchedCats


# Globals
class PlotKwargs(object):
    def __init__(self):
        self.ax= dict(fontweight='bold',fontsize='medium')
        self.text=dict(fontweight='bold',fontsize='medium',va='top',ha='left')
        self.leg=dict(frameon=True,fontsize='x-small')
        self.save= dict(bbox_inches='tight',dpi=150)

kwargs= PlotKwargs()

# Helpful funcs
def plot_all(obj, *args, **kwargs):
    for name in dir(obj):
        if name.startswith('plot'):
            attr = getattr(obj, name)
            attr(*args, **kwargs)


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

# Kaylans plots
class Kaylans(object):
    def __init__(self,ref_cat,obs_cat,imatch,\
                 ref_name='ref',obs_name='obs',savefig=False,\
                 plot_all=True):
        ##plot_all(self,ref_cat,obs_cat,\
        ##plot_all(self,ref_cat,obs_cat,\
        ##        ref_name='ref',obs_name='obs')
        if plot_all:
            self.nobs(ref_cat,savefig=savefig) 
            self.sn_vs_mag(ref_cat, mag_minmax=(18.,26.),savefig=savefig)
            self.barchart(ref_cat[ imatch['ref'] ],obs_cat[ imatch['obs'] ],\
                          ref_name=ref_name,obs_name=obs_name,savefig=savefig,prefix='alldata')
            self.radec(ref_cat[ imatch['ref'] ],obs_cat[ imatch['obs'] ],\
                       ref_name=ref_name,obs_name=obs_name,savefig=savefig)
            self.confusion_matrix(ref_cat[ imatch['ref'] ],obs_cat[ imatch['obs'] ],\
                                   ref_name=ref_name,obs_name=obs_name,savefig=savefig)
            for mytype in ['PSF','SIMP','DEV','COMP']:
                self.stacked_confusion_matrix(ref_cat[ imatch['ref'] ],obs_cat[ imatch['obs'] ],\
                                       ref_name=ref_name,obs_name=obs_name,savefig=savefig,\
                                       band='z',mytype=mytype)
            self.delta_mag_vs_mag(ref_cat[ imatch['ref'] ],obs_cat[ imatch['obs'] ],\
                                  ref_name=ref_name,obs_name=obs_name,savefig=savefig,ylim=[-0.2,0.2])
            self.chi_v_gaussian(ref_cat[ imatch['ref'] ],obs_cat[ imatch['obs'] ],\
                                 low=-4.,hi=4.,savefig=savefig) 
        
    
    def nobs(self,tractor, prefix='',savefig=False):
        '''make histograms of nobs so can compare depths of g,r,z between the two catalogues
        tractor -- Tractor catalogue in a table'''   
        hi= np.max(tractor.get('decam_nobs')[:,[1,2,4]])
        fig,ax= plt.subplots(3,1)
        for i, band,iband in zip(range(3),['g','r','z'],[1,2,4]):
            ax[i].hist(tractor.get('decam_nobs')[:,iband],\
                       bins=hi+1,normed=True,cumulative=True,align='mid')
            xlab=ax[i].set_xlabel('nobs %s' % band, **kwargs.ax)
            ylab=ax[i].set_ylabel('CDF', **kwargs.ax)
        if savefig == True:
            plt.savefig('%snobs.png' % prefix, bbox_extra_artists=[xlab,ylab], **kwargs.save)
            plt.close()

    def radec(self,ref_tractor,obs_tractor,ref_name='ref',obs_name='obs',\
              prefix='',savefig=False): 
        '''ra,dec distribution of objects
        obj -- Single_TractorCat()'''
        fig,ax= plt.subplots(1,2,figsize=(10,5))
        plt.subplots_adjust(wspace=0.1)
        for cnt,cat,name in zip(range(2),[ref_tractor,obs_tractor],[ref_name,obs_name]):
            ax[cnt].scatter(cat.get('ra'), cat.get('dec'), \
                            edgecolor='b',c='none',lw=1.)
            xlab=ax[cnt].set_xlabel('RA', **kwargs.ax)
            ti=ax[cnt].set_title(name, **kwargs.ax)
        ylab=ax[0].set_ylabel('DEC', **kwargs.ax)
        xmin,xmax= min(ax[0].get_xlim()[0],ax[1].get_xlim()[0]), max(ax[0].get_xlim()[1],ax[1].get_xlim()[1])
        for cnt in range(2):
            ax[cnt].set_xlim([xmin,xmax])
        if savefig == True:
            plt.savefig('%sradec.png' % prefix,bbox_extra_artists=[xlab,ylab], **kwargs.save)
            plt.close()

    def barchart(self,ref_tractor,obs_tractor, prefix='',\
                       ref_name='ref',obs_name='test',savefig=False):
        '''bar chart'''
        c1= 'b'
        c2= 'r'
        ###
        types= ['PSF','SIMP','EXP','DEV','COMP']
        ind = np.arange(len(types))  # the x locations for the groups
        width = 0.35       # the width of the bars
        ###
        ht_ref, ht_obs= np.zeros(5,dtype=int),np.zeros(5,dtype=int)
        for cnt,typ in enumerate(types):
            ht_ref[cnt]= np.where(ref_tractor.get('type') == typ)[0].shape[0]
            ht_obs[cnt]= np.where(obs_tractor.get('type') == typ)[0].shape[0]
        ###
        fig, ax = plt.subplots()
        rects1 = ax.bar(ind, ht_ref, width, color=c1)
        rects2 = ax.bar(ind + width, ht_obs, width, color=c2)
        ylab= ax.set_ylabel("N")
        ti= ax.set_title("Total: %s=%d, %s=%d" % (ref_name,len(ref_tractor),obs_name,len(obs_tractor)))
        ax.set_xticks(ind + width)
        ax.set_xticklabels(types)
        ax.legend((rects1[0], rects2[0]), (ref_name,obs_name))
        if savefig:
            plt.savefig('barchart%s.png' % prefix, bbox_extra_artists=[ylab,ti], bbox_inches='tight',dpi=200)
            plt.close()

    def sn_vs_mag(self,tractor, mag_minmax=(18.,26.),prefix='',savefig=False):
        '''plots Signal to Noise vs. mag for each band'''
        min,max= mag_minmax
        # Bin up SN values
        bin_SN={}
        for band,iband in zip(['g','r','z'],[1,2,4]):
            bin_SN[band]= bin_up(tractor.get('decam_mag_nodust')[:,iband], \
                           tractor.get('decam_flux')[:,iband]*np.sqrt(tractor.get('decam_flux_ivar')[:,iband]),\
                           bin_minmax=mag_minmax)
        #setup plot
        fig,ax=plt.subplots(1,3,figsize=(9,3),sharey=True)
        plt.subplots_adjust(wspace=0.1)
        #plot SN
        for cnt,band,color in zip(range(3),['g','r','z'],['g','r','m']):
            #horiz line at SN = 5
            ax[cnt].plot([mag_minmax[0],mag_minmax[1]],[5,5],'k--',lw=2,label='S/N=5')
            #ax[2].text(26,5,'S/N = 5  ',**kwargs.text)
            #data
            ax[cnt].plot(bin_SN[band]['binc'], bin_SN[band]['q50'],c=color,ls='-',lw=2,label='median')
            ax[cnt].fill_between(bin_SN[band]['binc'],bin_SN[band]['q25'],bin_SN[band]['q75'],\
                                 color=color,alpha=0.25,label='25/75 percentiles')
        #labels
        for cnt,band in zip(range(3),['g','r','z']):
            ax[cnt].set_yscale('log')
            xlab=ax[cnt].set_xlabel('%s' % band, **kwargs.ax)
            ax[cnt].set_ylim(1,100)
            ax[cnt].set_xlim(mag_minmax)
        ylab=ax[0].set_ylabel('S/N', **kwargs.ax)
        ax[2].legend(loc='lower right',fontsize='medium')
        if savefig == True:
            plt.savefig('%sSN.png' % prefix, \
                        bbox_extra_artists=[xlab,ylab], **kwargs.save)
            plt.close()

    def create_confusion_matrix(self,ref_tractor,test_tractor):
        '''compares MATCHED reference (truth) to test (prediction)
        ref_obj,test_obj -- reference,test Single_TractorCat()
        return 5x5 confusion matrix and colum/row names'''
        cm=np.zeros((5,5))-1
        types=['PSF','SIMP','EXP','DEV','COMP']
        for i_ref,ref_type in enumerate(types):
            cut_ref= np.where(ref_tractor.get('type') == ref_type)[0]
            #n_ref= ref_obj.number_not_masked(['current',ref_type.lower()])
            for i_test,test_type in enumerate(types):
                n_test= np.where(test_tractor.get('type')[ cut_ref ] == test_type)[0].size
                if cut_ref.size > 0: cm[i_ref,i_test]= float(n_test)/cut_ref.size
                else: cm[i_ref,i_test]= np.nan
        return cm,types


    def confusion_matrix(self,ref_tractor,test_tractor, prefix='',savefig=False,\
                         ref_name='ref',obs_name='test'):
        '''plot confusion matrix
        ref_obj,test_obj -- reference,test Single_TractorCat()'''
        cm,ticknames= self.create_confusion_matrix(ref_tractor,test_tractor)
        plt.imshow(cm, interpolation='nearest', cmap=plt.cm.Blues,vmin=0,vmax=1)
        cbar=plt.colorbar()
        plt.xticks(range(len(ticknames)), ticknames)
        plt.yticks(range(len(ticknames)), ticknames)
        ylab=plt.ylabel('True (%s)' % ref_name, **kwargs.ax)
        xlab=plt.xlabel('Predicted (%s)' % obs_name, **kwargs.ax)
        for row in range(len(ticknames)):
            for col in range(len(ticknames)):
                if np.isnan(cm[row,col]):
                    plt.text(col,row,'n/a',va='center',ha='center')
                elif cm[row,col] > 0.5:
                    plt.text(col,row,'%.2f' % cm[row,col],va='center',ha='center',color='yellow')
                else:
                    plt.text(col,row,'%.2f' % cm[row,col],va='center',ha='center',color='black')
        if savefig == True:
            plt.savefig('%sconfus_matrix.png' % prefix,bbox_extra_artists=[xlab,ylab], **kwargs.save)
            plt.close()

    def create_stack(self,ans=np.array(['PSF','PSF','PSF']),pred=np.array(['PSF','SIMP','COMP']),\
                     types=['PSF','SIMP','EXP','DEV','COMP'],slim=True):
        '''compares classifications of matched objects, returns 2D array which is conf matrix and xylabels
        return 5x5 confusion matrix and colum/row names
       answer_type,predict_type -- arrays of same length with reference and prediction types'''
        for typ in set(ans): assert(typ in types)
        for typ in set(pred): assert(typ in types)
        # if a type was not in answer (training) list then don't put in cm
        if slim: ans_types= set(ans)
        # put in cm regardless
        else: ans_types= set(types)
        cm=np.zeros((len(ans_types),len(types)))-1
        for i_ans,ans_type in enumerate(ans_types):
            ind= np.where(ans == ans_type)[0]
            for i_pred,pred_type in enumerate(types):
                n_pred= np.where(pred[ind] == pred_type)[0].size
                if ind.size > 0: cm[i_ans,i_pred]= float(n_pred)/ind.size # ind.size is constant for loop over pred_types
                else: cm[i_ans,i_pred]= np.nan
        if slim: return cm,ans_types,types #size ans_types != types
        else: return cm,types

    def plot_stack(self,cm_stack,stack_names,all_names, \
                   ref_name='ref',obs_name='test',\
                   prefix='',savefig=False):
        '''cm_stack -- list of single row confusion matrices
        stack_names -- list of same len as cm_stack, names for each row of cm_stack'''
        # combine list into single cm
        cm=np.zeros((len(cm_stack),len(all_names)))+np.nan
        for i in range(cm.shape[0]): cm[i,:]= cm_stack[i]
        # make usual cm, but labels repositioned
        plt.imshow(cm, interpolation='nearest', cmap=plt.cm.Blues)
        cbar=plt.colorbar()
        plt.xticks(range(len(all_names)), all_names)
        plt.yticks(range(len(stack_names)), stack_names)
        ylab=plt.ylabel('True (%s)' % ref_name)
        xlab=plt.xlabel('Predicted (%s)' % obs_name)
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
        if savefig == True:
            plt.savefig('%sconfus_matrix_stack.png' % prefix, bbox_extra_artists=[xlab,ylab], **kwargs.save)
            plt.close()

    def stacked_confusion_matrix(self,ref_matched,obs_matched,\
                                 ref_name='ref',obs_name='test', \
                                 outname='test.png',savefig=False,\
                                 band='z',mytype='PSF',prefix=''):
        cm_stack,stack_names=[],[]
        iband= dict(g=1,r=2,z=4)[band]
        for magmin,magmax in zip([18,20,22,23.5],[20,22,23.5,24.5]):
            band_and_type= (ref_matched.get('decam_mag_nodust')[:,iband] > magmin)* \
                           (ref_matched.get('decam_mag_nodust')[:,iband] <= magmax)*\
                           (ref_matched.get('type') == mytype)
            stack_names+= ["%s: %.1f < %s < %.1f" % (mytype,magmin,band,magmax)]
            cm,ans_names,all_names= self.create_stack(ans=ref_matched.get('type')[band_and_type],\
                                                      pred=obs_matched.get('type')[band_and_type])
            # If band_and_type is empty, fill with -1
            if len(cm) == 0: 
                cm,ans_names= np.zeros((1,5))-1, set([mytype])
            cm_stack+= [cm]
        self.plot_stack(cm_stack, stack_names,all_names, \
                   ref_name=ref_name,obs_name=obs_name,prefix='%s%s-%s-' % (prefix,mytype,band),savefig=savefig)


    def matched_dist(self,obj,dist, prefix=''):
        '''dist -- array of distances in degress between matched objects'''
        pixscale=dict(decam=0.25,bass=0.45)
        # Plot
        fig,ax=plt.subplots()
        ax.hist(dist*3600,bins=50,color='b',align='mid')
        ax2 = ax.twiny()
        ax2.hist(dist*3600./pixscale['bass'],bins=50,color='g',align='mid',\
                    visible=False)
        xlab= ax.set_xlabel("arcsec")
        xlab= ax2.set_xlabel("pixels [BASS]")
        ylab= ax.set_ylabel("Matched")
        if savefig == True:
            plt.savefig('%sdistance.png' % prefix, \
                        bbox_extra_artists=[xlab,ylab], **kwargs.save)
            plt.close()

    def chi_v_gaussian(self,ref_tractor,test_tractor,\
                       low=-8.,hi=8., prefix='',savefig=False):
        # Compute Chi
        chi={} 
        for band,iband in zip(['g','r','z'],[1,2,4]):
            chi[band]= (ref_tractor.get('decam_flux')[:,iband]-test_tractor.get('decam_flux')[:,iband])/\
                       np.sqrt( np.power(ref_tractor.get('decam_flux_ivar')[:,iband],-1)+\
                                np.power(test_tractor.get('decam_flux_ivar')[:,iband],-1))
        for b_low,b_hi in zip([21,22,23.5,20,18],[22,23.5,24.5,25,20]):
            #loop over mag bins, one 3 panel for each mag bin
            hist= dict(g=0,r=0,z=0)
            binc= dict(g=0,r=0,z=0)
            stats=dict(g=0,r=0,z=0)
            # Counts per bin
            for band,iband in zip(['g','r','z'],[1,2,4]):
                imag= np.all((b_low <= ref_tractor.get('decam_mag_nodust')[:,iband],\
                              ref_tractor.get('decam_mag_nodust')[:,iband] < b_hi),axis=0)
                hist[band],bins= np.histogram(chi[band][imag],\
                                        range=(low,hi),bins=50,normed=True)
                db= (bins[1:]-bins[:-1])/2
                binc[band]= bins[:-1]+db
            # Unit gaussian N(0,1)
            G= scipy.stats.norm(0,1)
            xvals= np.linspace(low,hi)
            # Plot for each mag range
            fig,ax=plt.subplots(3,1,figsize=(3,8),sharex=True)
            plt.subplots_adjust(hspace=0)
            for cnt,band in zip(range(3),['g','r','z']):
                ax[cnt].step(binc[band],hist[band], where='mid',c='b',lw=2) #label="%.1f < mag < %.1f" % (b_low,b_hi))
                ax[cnt].plot(xvals,G.pdf(xvals),c='g') #,label=r'$N(0,1)$')
            #labels
            for cnt,band in zip(range(3),['g','r','z']):
                ax[cnt].text(0.9,0.9,band,va='center',ha='center',transform=ax[cnt].transAxes)
                ax[cnt].set_ylim(0,0.6)
                ax[cnt].set_xlim(low,hi)
                ylab=ax[cnt].set_ylabel('PDF', **kwargs.ax)
            xlab=ax[2].set_xlabel(r'$(F_{1}-F_{2})/\sqrt{\sigma^2_{1}+\sigma^2_{1}}$', **kwargs.ax)
            ti=ax[0].set_title("%.1f < mag < %.1f" % (b_low,b_hi))
            ax[2].legend(loc='upper left',fontsize='medium')
            # Need unique name
            #name=os.path.basename(outname).replace('.ipynb','')+'_%d-%d.png' % (b_low,b_hi) 
            if savefig == True:
                plt.savefig('%schi_%d_%d.png' % (prefix,b_low,b_hi), \
                            bbox_extra_artists=[ti,xlab,ylab], **kwargs.save)
                plt.close()

    def delta_mag_vs_mag(self,ref_tractor,test_tractor, ref_name='ref',obs_name='test',\
                         prefix='',savefig=False,ylim=[-0.1,0.1]):
        fig,ax=plt.subplots(1,3,figsize=(12,3),sharey=True)
        plt.subplots_adjust(wspace=0.2)
        for cnt,iband in zip(range(3),[1,2,4]):
            delta= ref_tractor.get('decam_mag_nodust')[:,iband] - test_tractor.get('decam_mag_nodust')[:,iband]
            ax[cnt].scatter(ref_tractor.get('decam_mag_nodust')[:,iband],delta,\
                            c='none',facecolors='none',edgecolors='b',s=2,marker='o') 
        for cnt,band,maglim in zip(range(3),['g','r','z'],[24.5,23.5,22.]):
            xlab=ax[cnt].set_xlabel('%s (%s)' % (band,ref_name), **kwargs.ax)
            ax[cnt].set_ylim(ylim)
            ax[cnt].set_xlim(maglim-3,maglim+0.5)
        ax[cnt].yaxis.tick_right()
        ylab=ax[0].set_ylabel('mag(%s) - mag(%s)' % (ref_name,obs_name), **kwargs.ax)
        if savefig == True:
            plt.savefig('%sdelta_mag.png' % prefix, bbox_extra_artists=[xlab,ylab], **kwargs.save)
            plt.close()


    #def n_per_deg2(obj,deg2=1., req_mags=[24.,23.4,22.5],name=''):
    #    '''compute number density in each bin for each band mag [18,requirement]
    #    deg2 -- square degrees spanned by sources in obj table
    #    req_mags -- image requirements grz<=24,23.4,22.5'''
    #    bin_nd={}
    #    for band,iband,req in zip(['g','r','z'],[1,2,4],req_mags):
    #        bin_nd[band]={}
    #        bins= np.linspace(18.,req,num=15)
    #        bin_nd[band]['cnt'],junk= np.histogram(obj.t['decam_mag_nodust'][:,iband], bins=bins)
    #        bin_nd[band]['binc']= (bins[1:]+bins[:-1])/2.
    #        # bins and junk should be identical arrays
    #        assert( np.all(np.all((bins,junk),axis=0)) )
    #    # Plot
    #    fig,ax=plt.subplots(1,3,figsize=(9,3),sharey=True)
    #    plt.subplots_adjust(wspace=0.25)
    #    for cnt,band in zip(range(3),['g','r','z']):
    #        ax[cnt].step(bin_nd[band]['binc'],bin_nd[band]['cnt']/deg2, \
    #                     where='mid',c='b',lw=2)
    #    #labels
    #    for cnt,band in zip(range(3),['g','r','z']):
    #        xlab=ax[cnt].set_xlabel('%s' % band, **kwargs.ax) 
    #    ylab=ax[0].set_ylabel('counts/deg2', **kwargs.ax)
    #    # Make space for and rotate the x-axis tick labels
    #    fig.autofmt_xdate()
    #    plt.savefig(os.path.join(obj.outdir,'n_per_deg2_%s.png' % name), \
    #                bbox_extra_artists=[xlab,ylab], **kwargs.save)
    #    plt.close()


def noJunk(cat):
    '''cat is a fits_table object
    keep only S/N>5 '''
    i = (cat.get('brick_primary')) & (cat.get('decam_anymask')[:,1]==0) & \
        (cat.get('decam_anymask')[:,2]==0) & (cat.get('decam_anymask')[:,4]==0) & \
        (cat.get('decam_flux')[:,1]>0)  &  (cat.get('decam_flux')[:,2]>0) &  \
        (cat.get('decam_flux')[:,4]>0) & \
        (cat.get('decam_flux_ivar')[:,1]*(cat.get('decam_flux')[:,1])**2>25)  & \
        (cat.get('decam_flux_ivar')[:,2]*(cat.get('decam_flux')[:,2])**2>25)  & \
        (cat.get('decam_flux_ivar')[:,4]*(cat.get('decam_flux')[:,4])**2>25)
    return i

def areaCat(cat):
    '''cat is a fits_table object'''
    area = ( np.max(cat.get('ra')) - np.min(cat.get('ra')) ) * \
           ( np.max(cat.get('dec')) - np.min(cat.get('dec')) )*\
           np.cos( np.pi * np.mean(cat.get('dec')) / 180. )
    return area


class EnriqueCosmos(object):
    def __init__(self,ref_cat,obs_cat,imatch,\
                 ref_name='ref',obs_name='obs',\
                 savefig=False):
        self.compareMags(ref_cat[ imatch['ref'] ],obs_cat[ imatch['obs'] ],\
                         N1=ref_name,N2=obs_name,savefig=savefig)
    
    def compareMags(self,dr3,dr2,\
                   N1='ref',N2='obs',savefig=False):
        '''Enriques cosmso script used cosmo40/ as reference and cosmos41 as observed
        so N1=40 and N2=41
        catalague matching to the ref,
        he also called these dr3 for cosmos40, dr2 for cosmos41'''
        gfun = lambda x, m0, s0 : scipy.stats.norm.pdf(x,loc=m0,scale=s0)
        f=open("depth-comparison-new.txt","w")
        
        i_junk= noJunk(dr3)
        i_junk*= noJunk(dr2)
        dr2= dr2[i_junk]
        dr3= dr3[i_junk]
        areaDR3 = areaCat(dr3)
        areaDR2 = areaCat(dr2)

        #maglimits = dict(g=(21., 24.), r=(21., 23.5), z=(21., 22.5))
        # full-depth
        #maglimits = dict(g=(21.5, 24.), r=(21.5, 24.0), z=(21.5, 23))
        maglimits = dict(g=(23, 24.), r=(22.4, 23.4), z=(21.5, 22.5))

        # Check whether fitting the histogram introduces bias (binning is sinning)...
        # not much.
        G = np.random.normal(size=10000)
        nnn,bbb, ppp=plt.hist(G, bins=np.arange(-6,6,0.2), histtype='step', normed=True)
        out = curve_fit(gfun,(bbb[1:]+bbb[:-1])/2.,nnn,p0=(0,1))
        print('Fit to Gaussian random draws:', out)
        plt.clf()

        f.write(N1 + " " + N2 + '\n')

        bands= {"g","r","z"}

        for band in bands:
            wb ={"g":1,"r":2,"z":4}
            w = int(wb[band])
            print "band,w=", band,w,N1,N2
            mag_dr2 = dr2.get('decam_mag_nodust')
            mag_dr3 = dr3.get('decam_mag_nodust')
            mag_mean = 22.5 - 2.5 * np.log10((dr2.get('decam_flux')[:,w] + dr3.get('decam_flux')[:,w])/2.)

            iv_dr2 = dr2.get('decam_flux_ivar')[:,w]
            iv_dr3 = dr3.get('decam_flux_ivar')[:,w]

            df = dr3.get('decam_flux')[:,w] - dr2.get('decam_flux')[:,w]
            sigma = (1./dr2.get('decam_flux_ivar')[:,w] + 1./dr3.get('decam_flux_ivar')[:,w])**(0.5)

            plt.figure(2,(5,5))
            plt.axes([0.17,0.15,0.75,0.75])
            plt.plot(np.arange(-6,6,0.2), scipy.stats.norm.pdf(np.arange(-6,6,0.2),loc=0,scale=1.), 'k--', lw=0.5, label='N(0,1)')
            maglo,maghi = maglimits[band]

            ok = ((mag_mean > maglo) & (mag_mean < maghi) &
            (iv_dr2 > 0) & (iv_dr3 > 0))

            wb ={"g":1,"r":2,"z":4}
            w = int(wb[band])
            xwb= wb
            del xwb[band]
            xbands=xwb.keys()
            for xband in xbands:   # apply magnitude cuts in all bands
                wx= int(xwb[xband])
                mag_mean2 = 22.5 - 2.5 * np.log10((dr2.get('decam_flux')[:,wx] + \
                                                   dr3.get('decam_flux')[:,wx])/2.)
                maglo2,maghi2 = maglimits[xband]
                ok = (ok) & (mag_mean2 > maglo2) & (mag_mean2 < maghi2)
            
            sig = df[ok] / sigma[ok]
            mean_df=sum(sig)/len(sig)
            sigma_df=(sum((sig-mean_df)**2)/len(sig))**(0.5)

            xcor = np.zeros(len(xwb))
            ib= 0
            for xband in xbands:  # estimate covariance in xcor
                wx= int(xwb[xband])          
                dfx = dr3.get('decam_flux')[:,wx] - dr2.get('decam_flux')[:,wx]
                sigmax = (1./dr2.get('decam_flux_ivar')[:,wx] + 1./dr3.get('decam_flux_ivar')[:,wx])**(0.5)
                sigx = dfx[ok] / sigmax[ok]
                meanx_df=sum(sigx)/len(sigx)
                sigmax_df=(sum((sigx-meanx_df)**2)/len(sigx))**(0.5)
                xcor[ib]=sum((sig-mean_df)*(sigx-meanx_df))/len(sigx)/sigmax_df/sigma_df
                ib=ib+1

            xbandst=" Cov["+band+"-"+str(xbands[0])+","+band+"-"+str(xbands[1])+"]"
            txtr1=xbandst+"=["+str(np.round(xcor[0],2))+","+str(np.round(xcor[1],2))+"]"
            print txtr1
            out1= (mean_df,sigma_df)
            out3 = (np.median(sig), (np.percentile(sig,84) - np.percentile(sig,16))/2.)
            txt1="$\sigma$="+str(np.round(out1[1],2))+" , "+str(np.round(out3[1],2))

            nnn,bbb, ppp=plt.hist(sig, bins=np.arange(-6,6,0.2),histtype='step',label="All "+txt1+" "+txtr1, normed=True)
            out = curve_fit(gfun,(bbb[1:]+bbb[:-1])/2.,nnn,p0=(0,1))
            print  str(out[0]),txt1
            #gfit=scipy.stats.norm.pdf((bbb[1:]+bbb[:-1])/2.,loc=out[0][0],scale=out[0][1])
            #gfit1=gfit[(abs(bbb)<2)]
            #nnn1=nnn[(abs(bbb)<2)]
            #chi2=sum((nnn1-gfit1)**2/nnn1)*len(sig) # error is sqrt(N)/len(sig)
            #print "chi2=",chi2,chi2/(len(nnn1)-3)
            kst=sm.stats.lillifors(sig) # https://en.wikipedia.org/wiki/Lilliefors_test

            plt.plot(np.arange(-6,6,0.2), scipy.stats.norm.pdf(np.arange(-6,6,0.2),loc=out[0][0],scale=out[0][1]), 'b--', lw=2, label='All N fit [mean,$\sigma$]='+str(np.round(out[0],2))+' logP='+str(np.round(np.log10(kst[1]),1)))

            #ok2 = (ok) & (dr2.get('type')=="PSF") & (dr3.get('type')=='PSF')
            ok2 = (ok) & (dr2.get('type')=="SIMP") & (dr3.get('type')=='SIMP')

            sig2 = df[ok2] / sigma[ok2]
            mean_df2=sum(sig2)/len(sig2)
            sigma_df2=(sum((sig2-mean_df2)**2)/len(sig2))**(0.5)
            #sigma2_df=(np.var(sig))**(0.5)
            ib= 0
            for xband in xbands:
                wx= int(xwb[xband])
                dfx = dr3.get('decam_flux')[:,wx] - dr2.get('decam_flux')[:,wx]
                sigmax = (1./dr2.get('decam_flux_ivar')[:,wx] + 1./dr3.get('decam_flux_ivar')[:,wx])**(0.5)
                sigx = dfx[ok2] / sigmax[ok2]
                meanx_df=sum(sigx)/len(sigx)
                sigmax_df=(sum((sigx-meanx_df)**2)/len(sigx))**(0.5)
                xcor[ib]=sum((sig2-mean_df2)*(sigx-meanx_df))/len(sigx)/sigmax_df/sigma_df2
                ib=ib+1
            txtr2=xbandst+"=["+str(np.round(xcor[0],2))+","+str(np.round(xcor[1],2))+"]"
            print txtr2

            out4= (mean_df2,sigma_df2)
            out5 = (np.median(sig2), (np.percentile(sig2,84) - np.percentile(sig2,16))/2.)
            #txt2="$\sigma$="+str(np.round(out4[1],2))+" s68="+str(np.round(out5[1],2))
            txt2="$\sigma$="+str(np.round(out4[1],2))+" , "+str(np.round(out5[1],2))
            
            nnn,bbb,ppp = plt.hist(sig2, bins=np.arange(-6,6,0.2), histtype='step', label='SIMP '+txt2+" "+txtr2, normed=True)
            out2 = curve_fit(gfun,(bbb[1:]+bbb[:-1])/2.,nnn,p0=(0,1))
            print str(out2[0]),txt2

            #gfit=scipy.stats.norm.pdf((bbb[1:]+bbb[:-1])/2.,loc=out2[0][0],scale=out2[0][1])
            #gfit2=gfit[(abs(bbb)<2)]
            #nnn2=nnn[(abs(bbb)<2)]
            #chi22=sum((nnn2-gfit2)**2/nnn2)*len(sig2) # error is 1/sqrt(N)
            #print "chi22=",chi22,chi22/(len(nnn2)-3)     
            kst2=sm.stats.lillifors(sig2) #https://en.wikipedia.org/wiki/Lilliefors_test

            plt.plot(np.arange(-6,6,0.2), scipy.stats.norm.pdf(np.arange(-6,6,0.2),\
                     loc=out2[0][0],scale=out2[0][1]), 'g--', lw=2, \
                     label=r'SIMP N fit [mean,$\sigma$]='+\
                           str(np.round(out2[0],2))+\
                           ' logP='+str(np.round(np.log10(kst2[1]),1)))

            print  str(out2[0]),out4,out5

            f.write(band + str(out[0])+str(out1)+str(out3)+str(out2[0])+str(out4)+str(out5))
            f.write('\n')

            #ok = (dr2.get('type')=="PSF")&(dr3['decam_psfsize'].T[1]<1.5)
            #plt.hist(df[ok] / sigma[ok], bins=np.arange(-6,6,0.1), weights = np.ones_like(df[ok])/areaDR3, histtype='step', label='type PSF, PSFsize<1.5', normed=True)
            plt.xlabel("[ "+band+"("+N1+")-"+band+"("+N2+") ] / sqrt[ Var("+N2+")+Var("+N1+") ]")
            plt.ylabel('Normed counts %g < %s < %g' %(maglo,band,maghi))
            plt.xlim((-4,4))
            plt.ylim((0,0.7))
            gp = plt.legend(loc=2, fontsize=10)
            gp.set_frame_on(False)
            plt.title('Cosmos '+band+N1+"-"+N2+" #= "+str(len(sig))+'(All) '+str(len(sig2)))
            plt.grid()
            if savefig == True:
                plt.savefig(N1+'-'+N2+"-enrique-cosmos-"+band+".png")
                plt.clf()

class Dustins(object):
    def __init__(self,ref_cat,obs_cat,imatch,d2d,\
                 ref_name='ref',obs_name='obs',savefig=False,\
                 plot_all=True):
        #plot_all(self,ref_cat,obs_cat,
        #        ref_name='ref',obs_name='obs')
        self.outdir='./'
        self.cuts= Cuts4MatchedCats(ref_cat[ imatch['ref']], obs_cat[imatch['obs']])
        # Plot
        if plot_all:
            self.match_distance(d2d,savefig=savefig)
            self.fluxflux(ref_cat[ imatch['ref'] ],obs_cat[ imatch['obs'] ],\
                                      ref_name=ref_name,obs_name=obs_name,savefig=savefig)
            self.fluxdeltaflux(ref_cat[ imatch['ref'] ],obs_cat[ imatch['obs'] ],\
                                      ref_name=ref_name,obs_name=obs_name,savefig=savefig)
            self.grz_ratio_fluxerr(ref_cat[ imatch['ref'] ],obs_cat[ imatch['obs'] ],\
                                      ref_name=ref_name,obs_name=obs_name,savefig=savefig)
            self.magmag(ref_cat[ imatch['ref'] ],obs_cat[ imatch['obs'] ],\
                                      ref_name=ref_name,obs_name=obs_name,savefig=savefig)
            self.magdeltamag(ref_cat[ imatch['ref'] ],obs_cat[ imatch['obs'] ],\
                                      ref_name=ref_name,obs_name=obs_name,savefig=savefig)
            self.deltamag_err(ref_cat[ imatch['ref'] ],obs_cat[ imatch['obs'] ],\
                                      ref_name=ref_name,obs_name=obs_name,savefig=savefig)
            self.deltamag_errbars(ref_cat[ imatch['ref'] ],obs_cat[ imatch['obs'] ],\
                                      ref_name=ref_name,obs_name=obs_name,savefig=savefig)
            self.stephist(ref_cat[ imatch['ref'] ],obs_cat[ imatch['obs'] ],\
                                      ref_name=ref_name,obs_name=obs_name,savefig=savefig)
            #self.curvehist(ref_cat[ imatch['ref'] ],obs_cat[ imatch['obs'] ],\
            #                          ref_name=ref_name,obs_name=obs_name,savefig=savefig)

    def match_distance(self,dist,range=(0,1),prefix='',savefig=False):
        fig,ax=plt.subplots()
        ax.hist(dist * 3600., 100,range=range)
        third_pix= (1./3)*0.262
        ax.plot([third_pix]*2, ax.get_ylim(),'k--')
        ax.text(1.05*third_pix, 0.9*ax.get_ylim()[1],'1/3 pixel',ha='left',fontsize='medium')
        plt.xlabel('Match distance (arcsec)')
        if savefig == True:
            plt.savefig(os.path.join(self.outdir,'%smatch_dist.png' % prefix),bbox_inches='tight',dpi=150)
            plt.close()
    
    def fluxflux(self,matched1,matched2,\
                 ref_name='ref',obs_name='obs',prefix='',savefig=False):
        fig,ax= plt.subplots(1,3)
        for cnt,iband,band,cc in [(0,1,'g','g'),(1,2,'r','r'),(2,4,'z','m')]:
            K = np.flatnonzero(self.cuts.good[band]) 

            print('Median mw_trans', band, 'is',
                  np.median(matched1.decam_mw_transmission[:,iband]))
            ax[cnt].errorbar(matched1.decam_flux[K,iband],
                         matched2.decam_flux[K,iband],
                         fmt='.', color=cc,
                         xerr=1./np.sqrt(matched1.decam_flux_ivar[K,iband]),
                         yerr=1./np.sqrt(matched2.decam_flux_ivar[K,iband]),
                         alpha=0.1,
                         )
            ax[cnt].set_xlabel('%s flux: %s' % (ref_name, band))
            ax[cnt].set_ylabel('%s flux: %s' % (obs_name, band))
            ax[cnt].plot([-1e6, 1e6], [-1e6,1e6], 'k-', alpha=1.)
            ax[cnt].axis([-100, 1000, -100, 1000])
        if savefig == True:
            plt.savefig(os.path.join(self.outdir,'%sfluxflux_%s.png' % (prefix,band)))
            plt.close()

    def fluxdeltaflux(self,matched1,matched2,\
                     ref_name='ref',obs_name='obs',prefix='',savefig=False):
        fig,ax= plt.subplots(1,3)
        for cnt,iband,band,cc in [(0,1,'g','g'),(1,2,'r','r'),(2,4,'z','m')]:
            K = np.flatnonzero(self.cuts.good[band])
            P = np.flatnonzero(self.cuts.good[band] * self.cuts.psf1 * self.cuts.psf2)  

            mag1, magerr1 = matched1.decam_mag_nodust[:,iband],1./np.sqrt(matched1.decam_mag_ivar_nodust[:,iband])

            iv1 = matched1.decam_flux_ivar[:, iband]
            iv2 = matched2.decam_flux_ivar[:, iband]
            std = np.sqrt(1./iv1 + 1./iv2)

            ax[cnt].plot(mag1[K],
                     (matched2.decam_flux[K,iband] - matched1.decam_flux[K,iband]) / std[K],
                     '.', alpha=0.1, color=cc)
            ax[cnt].plot(mag1[P],
                     (matched2.decam_flux[P,iband] - matched1.decam_flux[P,iband]) / std[P],
                     '.', alpha=0.1, color='k')
            ax[cnt].set_ylabel('(%s - %s) flux / flux errors (sigma): %s' % \
                               (obs_name, ref_name, band))
            ax[cnt].set_xlabel('%s mag: %s' % (ref_name, band))
            ax[cnt].axhline(0, color='k', alpha=0.5)
            ax[cnt].axis([24, 16, -10, 10])
        if savefig == True:
            plt.savefig(os.path.join(self.outdir,'%sfluxdeltaflux_%s.png' % (prefix,band)))
            plt.close()

    def grz_ratio_fluxerr(self,matched1,matched2,\
                          ref_name='ref',obs_name='obs',prefix='',savefig=False):
        fig,ax= plt.subplots(1,3)
        lp,lt = [],[]
        for cnt,iband,band,cc in [(0,1,'g','g'),(1,2,'r','r'),(2,4,'z','m')]:        
            mag1, magerr1 = matched1.decam_mag_nodust[:,iband],1./np.sqrt(matched1.decam_mag_ivar_nodust[:,iband])

            iv1 = matched1.decam_flux_ivar[:, iband]
            iv2 = matched2.decam_flux_ivar[:, iband]
            std = np.sqrt(1./iv1 + 1./iv2)
            #std = np.hypot(std, 0.01)
            G= np.flatnonzero(self.cuts.good[band] * self.cuts.psf1 * self.cuts.psf2 *\
                              np.isfinite(mag1) *\
                              (mag1 >= 20) * (mag1 < dict(g=24, r=23.5, z=22.5)[band]))

            n,b,p = ax[cnt].hist((matched2.decam_flux[G,iband] -
                              matched1.decam_flux[G,iband]) / std[G],
                     range=(-4, 4), bins=50, histtype='step', color=cc,
                     normed=True)

            sig = (matched2.decam_flux[G,iband] -
                   matched1.decam_flux[G,iband]) / std[G]
            print('Raw mean and std of points:', np.mean(sig), np.std(sig))
            med = np.median(sig)
            print('sig= ',sig,'len(sig)=',len(sig))
            if len(sig) > 0:
                rsigma = (np.percentile(sig, 84) - np.percentile(sig, 16)) / 2.
            else: rsigma=-1
            print('Median and percentile-based sigma:', med, rsigma)
            lp.append(p[0])
            lt.append('%s: %.2f +- %.2f' % (band, med, rsigma))

        bins = []
        gaussint = []
        for blo,bhi in zip(b, b[1:]):
            c = scipy.stats.norm.cdf(bhi) - scipy.stats.norm.cdf(blo)
            c /= (bhi - blo)
            #bins.extend([blo,bhi])
            #gaussint.extend([c,c])
            bins.append((blo+bhi)/2.)
            gaussint.append(c)
        ax[cnt].plot(bins, gaussint, 'k-', lw=2, alpha=0.5)
        ax[cnt].set_xlabel('Flux difference / error (sigma)')
        ax[cnt].axvline(0, color='k', alpha=0.1)
        ax[cnt].set_ylim(0, 0.45)
        ax[cnt].legend(lp, lt, loc='upper right')
        if savefig == True:
            plt.savefig(os.path.join(self.outdir,'%sgrz_ratio_fluxerr.png' % prefix))
            plt.close()

    def magmag(self,matched1,matched2,\
              ref_name='ref',obs_name='obs',prefix='',savefig=False):
        fig,ax= plt.subplots(1,3)
        for cnt,iband,band,cc in [(0,1,'g','g'),(1,2,'r','r'),(2,4,'z','m')]:
            mag1, magerr1 = NanoMaggies.fluxErrorsToMagErrors(
                matched1.decam_flux[:,iband], matched1.decam_flux_ivar[:,iband])
            mag2, magerr2 = NanoMaggies.fluxErrorsToMagErrors(
                matched2.decam_flux[:,iband], matched2.decam_flux_ivar[:,iband])

            meanmag = NanoMaggies.nanomaggiesToMag((
                matched1.decam_flux[:,iband] + matched2.decam_flux[:,iband]) / 2.)

            K = np.flatnonzero(self.cuts.good[band])
            P = np.flatnonzero(self.cuts.good[band] * self.cuts.psf1 * self.cuts.psf2)

            ax[cnt].errorbar(mag1[K], mag2[K], fmt='.', color=cc,
                         xerr=magerr1[K], yerr=magerr2[K], alpha=0.1)
            ax[cnt].plot(mag1[P], mag2[P], 'k.', alpha=0.5)
            ax[cnt].set_xlabel('%s %s (mag)' % (ref_name, band))
            ax[cnt].set_ylabel('%s %s (mag)' % (obs_name, band))
            ax[cnt].plot([-1e6, 1e6], [-1e6,1e6], 'k-', alpha=1.)
            ax[cnt].axis([24, 16, 24, 16])
        if savefig == True:
            plt.savefig(os.path.join(self.outdir,'%smagmag_%s.png' % (prefix,band)))
            plt.close()

    def magdeltamag(self,matched1,matched2,\
              ref_name='ref',obs_name='obs',prefix='',savefig=False):
        fig,ax= plt.subplots(1,3)
        for cnt,iband,band,cc in [(0,1,'g','g'),(1,2,'r','r'),(2,4,'z','m')]:
            mag1, magerr1 = NanoMaggies.fluxErrorsToMagErrors(
                matched1.decam_flux[:,iband], matched1.decam_flux_ivar[:,iband])
            mag2, magerr2 = NanoMaggies.fluxErrorsToMagErrors(
                matched2.decam_flux[:,iband], matched2.decam_flux_ivar[:,iband])

            meanmag = NanoMaggies.nanomaggiesToMag((
                matched1.decam_flux[:,iband] + matched2.decam_flux[:,iband]) / 2.)

            K = np.flatnonzero(self.cuts.good[band])
            P = np.flatnonzero(self.cuts.good[band] * self.cuts.psf1 * self.cuts.psf2)

            ax[cnt].errorbar(mag1[K], mag2[K] - mag1[K], fmt='.', color=cc,
                         xerr=magerr1[K], yerr=magerr2[K], alpha=0.1)
            ax[cnt].plot(mag1[P], mag2[P] - mag1[P], 'k.', alpha=0.5)
            ax[cnt].set_xlabel('%s %s (mag)' % (ref_name, band))
            ax[cnt].set_ylabel('%s %s - %s %s (mag)' % (obs_name, band, ref_name, band))
            ax[cnt].axhline(0., color='k', alpha=1.)
            ax[cnt].axis([24, 16, -1, 1])
        if savefig == True:
            plt.savefig(os.path.join(self.outdir,'%smagdeltamag_%s.png' % (prefix,band)))
            plt.close()

    def deltamag_err(self,matched1,matched2,\
                     ref_name='ref',obs_name='obs',prefix='',savefig=False):
        fig,ax= plt.subplots(1,3)
        for cnt,iband,band,cc in [(0,1,'g','g'),(1,2,'r','r'),(2,4,'z','m')]:
            mag1, magerr1 = NanoMaggies.fluxErrorsToMagErrors(
                matched1.decam_flux[:,iband], matched1.decam_flux_ivar[:,iband])
            mag2, magerr2 = NanoMaggies.fluxErrorsToMagErrors(
                matched2.decam_flux[:,iband], matched2.decam_flux_ivar[:,iband])

            meanmag = NanoMaggies.nanomaggiesToMag((
                matched1.decam_flux[:,iband] + matched2.decam_flux[:,iband]) / 2.)

            K = np.flatnonzero(self.cuts.good[band])
            P = np.flatnonzero(self.cuts.good[band] * self.cuts.psf1 * self.cuts.psf2)

            magbins = np.arange(16, 24.001, 0.5)
            ax[cnt].plot(mag1[K], (mag2[K]-mag1[K]) / np.hypot(magerr1[K], magerr2[K]),
                         '.', color=cc, alpha=0.1)
            ax[cnt].plot(mag1[P], (mag2[P]-mag1[P]) / np.hypot(magerr1[P], magerr2[P]),
                         'k.', alpha=0.5)

            ax[cnt].set_xlabel('%s %s (mag)' % (ref_name, band))
            ax[cnt].set_ylabel('(%s %s - %s %s) / errors (sigma)' %
                       (obs_name, band, ref_name, band))
            ax[cnt].axhline(0., color='k', alpha=1.)
            ax[cnt].axis([24, 16, -10, 10])
        if savefig == True:
            plt.savefig(os.path.join(self.outdir,'%sdeltamag_err_%s.png' % (prefix,band)))
            plt.close()

    def deltamag_errbars(self,matched1,matched2,\
              ref_name='ref',obs_name='obs',prefix='',savefig=False):
        fig,ax= plt.subplots(1,3)
        for cnt,iband,band,cc in [(0,1,'g','g'),(1,2,'r','r'),(2,4,'z','m')]:
            mag1, magerr1 = NanoMaggies.fluxErrorsToMagErrors(
                matched1.decam_flux[:,iband], matched1.decam_flux_ivar[:,iband])
            mag2, magerr2 = NanoMaggies.fluxErrorsToMagErrors(
                matched2.decam_flux[:,iband], matched2.decam_flux_ivar[:,iband])

            meanmag = NanoMaggies.nanomaggiesToMag((
                matched1.decam_flux[:,iband] + matched2.decam_flux[:,iband]) / 2.)

            K = np.flatnonzero(self.cuts.good[band])
            P = np.flatnonzero(self.cuts.good[band] * self.cuts.psf1 * self.cuts.psf2)

            magbins = np.arange(16, 24.001, 0.5)
            y = (mag2 - mag1) / np.hypot(magerr1, magerr2)
            ax[cnt].plot(meanmag[P], y[P], 'k.', alpha=0.1)

            midmag = []
            vals = np.zeros((len(magbins)-1, 5))
            median_err1 = []

            iqd_gauss = scipy.stats.norm.ppf(0.75) - scipy.stats.norm.ppf(0.25)

            # FIXME -- should we do some stats after taking off the mean difference?

            for bini,(mlo,mhi) in enumerate(zip(magbins, magbins[1:])):
                I = P[(meanmag[P] >= mlo) * (meanmag[P] < mhi)]
                midmag.append((mlo+mhi)/2.)
                median_err1.append(np.median(magerr1[I]))
                if len(I) == 0:
                    continue
                # median and +- 1 sigma quantiles
                ybin = y[I]
                vals[bini,0] = np.percentile(ybin, 16)
                vals[bini,1] = np.median(ybin)
                vals[bini,2] = np.percentile(ybin, 84)
                # +- 2 sigma quantiles
                vals[bini,3] = np.percentile(ybin, 2.3)
                vals[bini,4] = np.percentile(ybin, 97.7)

                iqd = np.percentile(ybin, 75) - np.percentile(ybin, 25)

                print('Mag bin', midmag[-1], ': IQD is factor', iqd / iqd_gauss,
                      'vs expected for Gaussian;', len(ybin), 'points')

                # if iqd > iqd_gauss:
                #     # What error adding in quadrature would you need to make the IQD match?
                #     err = median_err1[-1]
                #     target_err = err * (iqd / iqd_gauss)
                #     sys_err = np.sqrt(target_err**2 - err**2)
                #     print('--> add systematic error', sys_err)

            # ~ Johan's cuts
            mlo = 21.
            mhi = dict(g=24., r=23.5, z=22.5)[band]
            I = P[(meanmag[P] >= mlo) * (meanmag[P] < mhi)]
            print('y=',y)
            print('I=',I)
            ybin = y[I]
            print('ybin =',ybin)
            if len(ybin) > 0:
                iqd = np.percentile(ybin, 75) - np.percentile(ybin, 25)
            else: iqd=-1
            print('Mag bin', mlo, mhi, 'band', band, ': IQD is factor',
                  iqd / iqd_gauss, 'vs expected for Gaussian;', len(ybin), 'points')
            if iqd > iqd_gauss:
                # What error adding in quadrature would you need to make
                # the IQD match?
                err = np.median(np.hypot(magerr1[I], magerr2[I]))
                print('Median error (hypot):', err)
                target_err = err * (iqd / iqd_gauss)
                print('Target:', target_err)
                sys_err = np.sqrt((target_err**2 - err**2) / 2.)
                print('--> add systematic error', sys_err)

                # check...
                err_sys = np.hypot(np.hypot(magerr1, sys_err),
                                   np.hypot(magerr2, sys_err))
                ysys = (mag2 - mag1) / err_sys
                ysys = ysys[I]
                print('Resulting median error:', np.median(err_sys[I]))
                iqd_sys = np.percentile(ysys, 75) - np.percentile(ysys, 25)
                print('--> IQD', iqd_sys / iqd_gauss, 'vs Gaussian')
                # Hmmm, this doesn't work... totally overshoots.


            ax[cnt].errorbar(midmag, vals[:,1], fmt='o', color='b',
                         yerr=(vals[:,1]-vals[:,0], vals[:,2]-vals[:,1]),
                         capthick=3, zorder=20)
            ax[cnt].errorbar(midmag, vals[:,1], fmt='o', color='b',
                         yerr=(vals[:,1]-vals[:,3], vals[:,4]-vals[:,1]),
                         capthick=2, zorder=20)
            ax[cnt].axhline( 1., color='b', alpha=0.2)
            ax[cnt].axhline(-1., color='b', alpha=0.2)
            ax[cnt].axhline( 2., color='b', alpha=0.2)
            ax[cnt].axhline(-2., color='b', alpha=0.2)

            for mag,err,y in zip(midmag, median_err1, vals[:,3]):
                if not np.isfinite(err):
                    continue
                if y < -6:
                    continue
                plt.text(mag, y-0.1, '%.3f' % err, va='top', ha='center', color='k',
                         fontsize=10)

            ax[cnt].set_xlabel('(%s + %s)/2 %s (mag), PSFs' % (ref_name, obs_name, band))
            ax[cnt].set_ylabel('(%s %s - %s %s) / errors (sigma)' %
                       (obs_name, band, ref_name, band))
            ax[cnt].axhline(0., color='k', alpha=1.)

            ax[cnt].axvline(21, color='k', alpha=0.3)
            ax[cnt].axvline(dict(g=24, r=23.5, z=22.5)[band], color='k', alpha=0.3)

            ax[cnt].axis([24.1, 16, -6, 6])
        if savefig == True:
            plt.savefig(os.path.join(self.outdir,'%sdeltamag_errbars_%s.png' % (prefix,band)))
            plt.close()

    def stephist(self,matched1,matched2,\
                 ref_name='ref',obs_name='obs',prefix='',savefig=False):
        fig,ax= plt.subplots(1,3)
        for cnt,iband,band,cc in [(0,1,'g','g'),(1,2,'r','r'),(2,4,'z','m')]:
            mag1, magerr1 = NanoMaggies.fluxErrorsToMagErrors(
                matched1.decam_flux[:,iband], matched1.decam_flux_ivar[:,iband])
            mag2, magerr2 = NanoMaggies.fluxErrorsToMagErrors(
                matched2.decam_flux[:,iband], matched2.decam_flux_ivar[:,iband])

            meanmag = NanoMaggies.nanomaggiesToMag((
                matched1.decam_flux[:,iband] + matched2.decam_flux[:,iband]) / 2.)

            K = np.flatnonzero(self.cuts.good[band])
            P = np.flatnonzero(self.cuts.good[band] * self.cuts.psf1 * self.cuts.psf2)

            #magbins = np.append([16, 18], np.arange(20, 24.001, 0.5))
            if band == 'g':
                magbins = [20, 24]
            elif band == 'r':
                magbins = [20, 23.5]
            elif band == 'z':
                magbins = [20, 22.5]

            slo,shi = -5,5
            ha = dict(bins=25, range=(slo,shi), histtype='step', normed=True)
            y = (mag2 - mag1) / np.hypot(magerr1, magerr2)
            midmag = []
            nn = []
            rgbs = []
            lt,lp = [],[]
            I_empty=True
            for bini,(mlo,mhi) in enumerate(zip(magbins, magbins[1:])):
                I = P[(mag1[P] >= mlo) * (mag1[P] < mhi)]
                if len(I) == 0:
                    continue
                I_empty=False
                ybin = y[I]
                rgb = [0.,0.,0.]
                rgb[0] = float(bini) / (len(magbins)-1)
                rgb[2] = 1. - rgb[0]
                n,b,p = plt.hist(ybin, color=rgb, **ha)
                lt.append('mag %g to %g' % (mlo,mhi))
                lp.append(p[0])
                midmag.append((mlo+mhi)/2.)
                nn.append(n)
                rgbs.append(rgb)

            bins = []
            gaussint = []
            if I_empty == False:
                for blo,bhi in zip(b, b[1:]):
                    #midbin.append((blo+bhi)/2.)
                    #gaussint.append(scipy.stats.norm.cdf(bhi) -
                    #                scipy.stats.norm.cdf(blo))
                    c = scipy.stats.norm.cdf(bhi) - scipy.stats.norm.cdf(blo)
                    c /= (bhi - blo)
                    bins.extend([blo,bhi])
                    gaussint.extend([c,c])
            ax[cnt].plot(bins, gaussint, 'k-', lw=2, alpha=0.5)

            ax[cnt].legend(lp, lt)
            ax[cnt].set_xlim(slo,shi)
        if savefig == True:
            plt.savefig(os.path.join(self.outdir,'%sstephist_%s.png' % (prefix,band)))
            plt.close()

#    def curvehist(self,matched1,matched2,\
#              ref_name='ref',obs_name='obs',prefix='',savefig=False):
#        fig,ax= plt.subplots(1,3)
#        for cnt,iband,band,cc in [(0,1,'g','g'),(1,2,'r','r'),(2,4,'z','m')]:
#            mag1, magerr1 = NanoMaggies.fluxErrorsToMagErrors(
#                matched1.decam_flux[:,iband], matched1.decam_flux_ivar[:,iband])
#            mag2, magerr2 = NanoMaggies.fluxErrorsToMagErrors(
#                matched2.decam_flux[:,iband], matched2.decam_flux_ivar[:,iband])
#
#            meanmag = NanoMaggies.nanomaggiesToMag((
#                matched1.decam_flux[:,iband] + matched2.decam_flux[:,iband]) / 2.)
#
#            K = np.flatnonzero(self.cuts.good[band])
#            P = np.flatnonzero(self.cuts.good[band] * self.cuts.psf1 * self.cuts.psf2)
#
#            #magbins = np.append([16, 18], np.arange(20, 24.001, 0.5))
#            if band == 'g':
#                magbins = [20, 24]
#            elif band == 'r':
#                magbins = [20, 23.5]
#            elif band == 'z':
#                magbins = [20, 22.5]
#
#            slo,shi = -5,5
#            ha = dict(bins=25, range=(slo,shi), histtype='step', normed=True)
#            y = (mag2 - mag1) / np.hypot(magerr1, magerr2)
#            midmag = []
#            nn = []
#            rgbs = []
#            lt,lp = [],[]
#            I_empty=True
#            for bini,(mlo,mhi) in enumerate(zip(magbins, magbins[1:])):
#                I = P[(mag1[P] >= mlo) * (mag1[P] < mhi)]
#                if len(I) == 0:
#                    continue
#                I_empty=False
#                ybin = y[I]
#                rgb = [0.,0.,0.]
#                rgb[0] = float(bini) / (len(magbins)-1)
#                rgb[2] = 1. - rgb[0]
#                n,b,p = plt.hist(ybin, color=rgb, **ha)
#                lt.append('mag %g to %g' % (mlo,mhi))
#                lp.append(p[0])
#                midmag.append((mlo+mhi)/2.)
#                nn.append(n)
#                rgbs.append(rgb)
#
#            if I_empty == False:
#                bincenters = b[:-1] + (b[1]-b[0])/2.
#                lp = []
#                for n,rgb,mlo,mhi in zip(nn, rgbs, magbins, magbins[1:]):
#                    p = ax[cnt].plot(bincenters, n, '-', color=rgb)
#                    lp.append(p[0])
#                ax[cnt].plot(bincenters, gaussint[::2], 'k-', alpha=0.5, lw=2)
#                ax[cnt].legend(lp, lt)
#                ax[cnt].set_xlim(slo,shi)
#        if savefig == True:
#            plt.savefig(os.path.join(self.outdir,'%scurvehist_%s.png' % (prefix,band)))
#            plt.close()
