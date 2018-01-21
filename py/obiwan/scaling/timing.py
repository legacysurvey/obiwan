import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import os
import pickle
import json
import re
import pandas as pd

STAGES=['tims', 'mask_junk', 'srcs', 
        'fitblobs', 'coadds', 'writecat']


def add_scatter(ax,x,y,c='b',m='o',lab='hello',s=80,drawln=False):
	ax.scatter(x,y, s=s, lw=2.,facecolors='none',edgecolors=c, marker=m,label=lab)
	if drawln: ax.plot(x,y, c=c,ls='-')

class Plots(object):
    def __init__(self,tm):
        self.tm= tm

    def tractor_profile_plots(tm,name='tmp.png'):
        fig,ax=plt.subplots()
        xvals= np.arange(tm['stage'].size)+1
        print(tm['parallel'])
        add_scatter(ax,xvals, tm['serial']/60., c='b',m='o',lab='serial',drawln=True)
        add_scatter(ax,xvals, tm['parallel']/60., c='g',m='o',lab='parallel',drawln=True)
        plt.legend(loc='lower right',scatterpoints=1)
        #add_scatter(ax,xvals, tm['total']/60., c='b',m='o',lab='total')
        ax.set_xticks(xvals)
        ax.set_xticklabels(tm['stage'],rotation=45, ha='right')
        ax.set_yscale('log')
        ax.set_ylim([1e-3,1e2])
        xlab=ax.set_ylabel('Wall Time (min)')
        ylab=ax.set_xlabel('Tractor Stage')
        plt.savefig(name, bbox_extra_artists=[xlab,ylab], bbox_inches='tight',dpi=150)
        plt.close()

    def plot_wall_node(d):
        name='wall_v_nodes.png'
        fig,ax=plt.subplots()
        xvals= np.arange(d['nodes'].size)+1
        add_scatter(ax,xvals, d['tims_mean']/60., c='b',m='o',lab='tims',drawln=True)
        add_scatter(ax,xvals, d['fit_mean']/60., c='g',m='o',lab='fit',drawln=True)
        add_scatter(ax,xvals, d['tot_mean']/60., c='k',m='o',lab='total',drawln=True)
        plt.legend(loc='lower right',scatterpoints=1)
        #add_scatter(ax,xvals, tm['total']/60., c='b',m='o',lab='total')
        ax.set_xticks(xvals)
        names= np.zeros(d['nodes'].size).astype(str)
        for i in range(names.size): 
            names[i]= '%d/%d' % (d['cores'][i],d['nodes'][i])
        ax.set_xticklabels(names,rotation=45, ha='right')
        #ax.set_yscale('log')
        #ax.set_ylim([1e-3,1e3])
        ylab=ax.set_ylabel('Wall Time (min)')
        xlab=ax.set_xlabel('Cores/Nodes')
        plt.savefig(name, bbox_extra_artists=[xlab,ylab], bbox_inches='tight',dpi=150)
        plt.close()

def params_of_run(logfile):
    with open(logfile,'r') as f:
        bigstring= f.read()
    
    def get_param(expr,bigstring):
        a=re.search(expr,bigstring)
        return (bigstring[slice(a.regs[0][0],a.regs[0][1])]
                .split('=')[1]
                .replace(',',''))
    d={}
    d['rsdir']= get_param(r'rowstart=[0-9]+,',bigstring)
    d['nobj']= get_param(r'nobj=[0-9]+,',bigstring)
    d['brick']= get_param(r"brick='[0-9]{4}[mp][0-9]{3}',",bigstring).replace("'",'')
    return d

def number_injected(logfile,nobj=None):
    with open(logfile,'r') as f:
        bigstring= f.read()
    d={}
    a= re.search(r'INFO:decals_sim:sources.*?flagged as nearby [0-9]+?',bigstring)
    n_skip= (bigstring[slice(a.regs[0][0],a.regs[0][1])]
             .split(' ')[-1])
    d['frac_injected']= (nobj-int(n_skip))/float(nobj)
    return d

def time_per_stage(logfile):
    """Returns dict of seconds spend in each stage"""
    with open(logfile,'r') as f:
        bigstring= f.read()
    # rsdir
    a=re.search(r'rowstart=[0-9]+,',bigstring)
    rsdir= (bigstring[slice(a.regs[0][0],a.regs[0][1])]
            .split('=')[1]
            .replace(',',''))
    rsdir= (bigstring[slice(a.regs[0][0],a.regs[0][1])]
            .split('=')[1]
            .replace(',',''))

    t={}
    for stage in STAGES: 
        a=re.search(r'Resources for stage %s(.*\n)*?Grand total Wall:.*\n' % stage,
                    bigstring)
        print('stage=%s, a=' % stage,a)
        lines= bigstring[slice(a.regs[0][0],a.regs[0][1])].split('\n')
        print('lines=',lines)
        lines= pd.Series(lines) 
        line= lines[lines.str.contains('Grand total Wall')].str.split(r'\s+')
        assert(line.size == 1)
        assert(line.str[-1].values[0] == 'sec')
        t[stage]=line.str[-2].values[0]
    return t

def write_header(savenm):
    with open(savenm,'w') as foo:
        text= 'nobj brick rsdir frac_injected'
        for stage in STAGES:
            text += ' %s' % stage
        foo.write(text+'\n')
    print('Wrote header %s' % savenm)

def write_measurements(d,savenm='test.txt'):
    with open(savenm,'a') as foo:
        text= '%s %s %s %.3f' % (d['nobj'],d['brick'],d['rsdir'],d['frac_injected'])
        for stage in STAGES:
            text += ' %s' % d[stage]
        foo.write(text+'\n')
    print('Appended measurements %s' % savenm)
        

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description="test")
    parser.add_argument("--logfiles",action="store",required=True,
                        help="list of logfiles for the scaling run, e.g. nobj=500,1000,1500 for many bricks and rsdirs")
    parser.add_argument("--savenm",action="store",help='text file name to write measurements to',required=True)
    args = parser.parse_args()

    # Extract
    if not os.path.exists(args.savenm):
        write_header(args.savenm)
        fns= np.loadtxt(args.logfiles,dtype=str)
        for fn in fns:
            d= {**params_of_run(fn),
                **time_per_stage(fn)
               }
            d= {**d,
                **number_injected(fn,nobj=int(d['nobj']))
                }
            write_measurements(d, args.savenm)
    
    # Plots
    df= pd.read_csv(args.savenm,sep=' ')

