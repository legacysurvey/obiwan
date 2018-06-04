import matplotlib
matplotlib.use('Agg') # display backend
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import obiwan.qa.plots_common as plots

plt.rcParams['axes.labelsize'] = 14
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12

desi_depth= dict(g=24.0,
                 r=23.4,
                 z=22.5)

def fraction_recovered_csv(which,fn='fraction_recovered_csv.png',
                           glim=(20,26),rlim=(20,26),zlim=(20,26)):
    assert(which == 'cosmos')
    fig,axes=plt.subplots(1,3,figsize=(14,4))
    plt.subplots_adjust(wspace=0.05)
    xlim= dict(g=glim,
               r=rlim,
               z=zlim)

    get_fn= lambda subset,band: os.path.join('subset%d' % subset,
                                             'fraction_recovered_%s.csv' % band)
    for ax,band in zip(axes,'grz'):
        for subset,color in zip([60,64,69],'bgr'):
            df= pd.read_csv(get_fn(subset,band))
            ax.step(df['bins'],df['fraction'],where='mid',c=color,
                    label='subset %d' % subset)
        ax.axhline(0.5,c='k',ls='--')
        ax.axvline(desi_depth[band],c='k',ls='--')
        xlab=ax.set_xlabel('%s (AB mag, Truth)' % band)
    for ax in axes:
        ax.set_ylim(0.4,1)
    for ax in axes[1:]:
        ax.set_yticklabels([])
    ylab=axes[0].set_ylabel('Fraction Recovered')
    axes[0].legend(loc='upper left')
    plt.savefig(fn,bbox_extra_artists=[xlab,ylab], bbox_inches='tight')
    plt.close()
    print('Wrote %s' % fn)

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description='DECaLS simulations.')
    parser.add_argument('--which', choices=['cosmos','eboss','desi'],required=True)
    args = parser.parse_args()

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

    fraction_recovered_csv(args.which,**kw_lims)

