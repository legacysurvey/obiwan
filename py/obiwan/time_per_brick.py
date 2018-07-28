import matplotlib.pyplot as plt
import numpy as np
from astrometry.util.fits import fits_table

def add_scatter(ax,x,y,c='b',m='o',lab='hello',s=80,drawln=False):
    ax.scatter(x,y, s=s, lw=2.,facecolors='none',edgecolors=c, marker=m,label=lab)
    if drawln:
        ax.plot(x,y, c=c,ls='-')

def plot(T,name='stages.png'):
    fig,ax=plt.subplots()
    xvals= np.arange(len(stages))+1
    for stage in stages:
        add_scatter(ax,xvals,T.get(stage)/60., c='b',m='o',lab=stage,drawln=True)
    plt.legend(loc='upper right',scatterpoints=1)
    ax.set_xticks(xvals)
    ax.set_xticklabels(stages,rotation=45, ha='right')
    ax.set_yscale('log')
    xlab=ax.set_ylabel('Wall Time (min)')
    ylab=ax.set_xlabel('Tractor Stage')
    plt.savefig(name, bbox_extra_artists=[xlab,ylab], bbox_inches='tight',dpi=150)
    print('Wrote %s' % name)
    plt.close()

if __name__ == '__main__':
    stages= ['tims','mask_junk', 'srcs', 'fitblobs', 'coadds', 'wise_forced', 'writecat']

    T=fits_table()
    for stage in stages:
        b,_,t,_= np.loadtxt('obiwan/wall_%s.txt' % stage,dtype=str,unpack=True)
        if not 'brick' in T.get_columns():
            T.set('brick',b)
        T.set(stage,t.astype(np.float32))
    print('%d bricks' % len(T))

    # Report logs missing timing statement
    keep=np.ones(len(T),bool)
    for stage in stages:
        keep*= (np.isfinite(T.get(stage)))
    print('%d bricks missing at least one time print' % np.where(keep == False)[0].size)
    T.cut(keep)
    print('%d bricks remain' % len(T))

    plot(T)
