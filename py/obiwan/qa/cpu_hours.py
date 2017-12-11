import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def cpu_hrs(nodes,h,m,s, mpp=False,cori=True):
    """
    Args:
        nodes: number of nodes
        h,m,s: hours,min,sec
    """
    cores_per_node=32
    if not cori:
        cores_per_node=24
    mpp_factor= 1.
    if mpp:
        mpp_factor= 2.5
        if not cori:
            mpp_factor= 2.
    
    return mpp_factor * cores_per_node * nodes * h + m/60. + s/3600.

def dobash(cmd):
    print('UNIX cmd: %s' % cmd)
    if os.system(cmd): raise ValueError

if __name__ == '__main__':
    #sacct -A desi --user=kaylanb --format=JobID,State,NNodes,Elapsed -S 11/17/17 -E 12/03/17|grep COMPLETED > my_sacct.txt

    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--sacct_fn', type=str, required=True, help='output of sacct dumped to textfil')
    args = parser.parse_args()

    a=pd.read_csv(args.sacct_fn,
                  delim_whitespace=True,header=None,
                  names=['slurm_id','status','num_nodes','time'])

    for i,name in zip([0,1,2],['multi_hr','min','sec']):
        a[name]= a['time'].str.split(':').str[i]

    print(a['multi_hr'].str.split('-').str.len().value_counts())
    hasExtra24= a['multi_hr'].str.split('-').str.len() > 1
    a['extra_hrs']= np.zeros(a.shape[0])
    a.loc[hasExtra24,'extra_hrs']= 24*a.loc[hasExtra24,'multi_hr'].str.split('-').str[0].astype(float)
    a['extra_hrs'].value_counts()

    a['hrs']= np.zeros(a.shape[0])
    # No 01-05, just 05
    a.loc[~hasExtra24,'hrs']= a.loc[~hasExtra24,'multi_hr'].astype(float)
    # When 01-05, just take 05
    a.loc[hasExtra24,'hrs']= a.loc[hasExtra24,'multi_hr'].str.split('-').str[1].astype(float)
    a.loc[:,'hrs'] +=  a['extra_hrs']

    for col in ['min','sec']:
        a.loc[:,col]= a.loc[:,col].astype(float)

    a['cpu_hrs']= cpu_hrs(a['num_nodes'],a['hrs'],a['min'],a['sec'])
    print('total cpu hours (M):',a['cpu_hrs'].sum()/1e6)
    print('total MPP hours (M):',a['cpu_hrs'].sum()/1e6*2.5)


    ## Plots
    fig,ax= plt.subplots(2,2,figsize=(6,6))
    names=['num_nodes','hrs','min','sec']
    i=0
    for row in range(2):
        for col in range(2):
            sns.distplot(a[names[i]],ax=ax[row,col])
            i+=1
    fn='nodes_hrs_min_sec.png'
    plt.savefig(fn,dpi=150)
    print('Wrote %s' % fn)
