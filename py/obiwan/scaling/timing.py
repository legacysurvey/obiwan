import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import os
from argparse import ArgumentParser

def bash_result(cmd):
    res= subprocess.check_output(cmd,\
                        stderr=subprocess.STDOUT,\
                        shell=True)
    return res.strip()

def parse_mem(lines,nstages):
    d={}
    for key in ['stage','Wall','CPU','VmPeak','VmSize','VmRSS','VmData','maxrss']:
        d[key]=[]
    for line in lines:
        if 'Resources for stage' in line: d['stage']+= [line.split()[3]]
        elif 'Wall:' in line: 
            line= line.split(',')
			# Avoid extra worker thread columns if exists
            line= line[:len(d.keys())-1]
            for li in line:
                li=li.split()
                key=li[0][:-1]
                num=float(li[1])
#                 unit=li[2]
#                 print "key=-%s-, d.keys=" % key,d.keys()
                d[key]+= [num]
    for key in d.keys(): d[key]= np.array(d[key])
    assert(len(d['stage']) == nstages)
    return d

class Timing(object):
    def __init__(self):
        pass

def parse_time(lines,nstages):
	d= Timing()
	for key in ['stage','serial','parallel','total','util']:
		setattr(d,key,[])
	for line in lines:
		if 'Resources for stage' in line: 
			d['stage']+= [line.split()[3]]
		elif 'serial' in line:
			d['serial']+= [float(line.split()[3])]
		elif 'parallel' in line:
			d['parallel']+= [float(line.split()[3])]
		elif 'total Wall' in line:
			d['total']+= [float(line.split()[3])]
		elif 'CPU util' in line:
			d['util']+= [float(line.split()[4])]
	for key in d.keys(): d[key]= np.array(d[key])
	assert(len(d['stage']) == nstages)
	if d['stage'].size < d['serial'].size:
		# Serial, paralel, etc contains Grand Total vals
		for key in d.keys():
			if key != 'stage': d[key]= d[key][:d['stage'].size]
	return d

def parse_tractor_profile(fn):
    fobj=open(fn,'r')
    lines=fobj.readlines()
    fobj.close()
    lines=np.char.strip(lines)
    # one more item the funcs need
    nstages= int( bash_result("grep 'Resources for stage' %s | wc -l" % fn) )
    if 'mem_' in fn:
        return parse_mem(lines,nstages)
    elif 'time_' in fn:
        return parse_time(lines,nstages)
    else: raise ValueError

def add_scatter(ax,x,y,c='b',m='o',lab='hello',s=80,drawln=False):
	ax.scatter(x,y, s=s, lw=2.,facecolors='none',edgecolors=c, marker=m,label=lab)
	if drawln: ax.plot(x,y, c=c,ls='-')

def tractor_profile_plots(mem,tm,nthreads=1):
	name='time_v_stage_threads%d.png' % nthreads
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
    fobj=open(fn)
    tm[key]=pickle.load(fobj)
    fobj.close()
    # Plot
    fig,ax=plt.subplots()
    for key,col in zip(tm.keys(),['b','g']):
        xvals= np.arange(tm[key]['stage'].size)+1
        add_scatter(ax,xvals, (tm[key]['serial']+tm[key]['parallel'])/60., c=col,m='o',lab=key,drawln=True)
        #add_scatter(ax,xvals, tm['total']/60., c='b',m='o',lab='total')
    plt.legend(loc='upper right',scatterpoints=1)
    ax.set_xticks(xvals)
    ax.set_xticklabels(tm[key]['stage'],rotation=45, ha='right')
    ax.set_yscale('log')
    #ax.set_ylim([1e-2,1e1])
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



if __name__ == '__main__':
    # Tractor stdout file, parse profiling info
    parser = ArgumentParser(description="test")
    parser.add_argument("--logfile",action="store",required=True)
    parser.add_argument("--rsdir",action="store",type=str,default=0,required=True)
    parser.add_argument("--outdir",action="store",help='where to write outputs',default='.',required=False)
    args = parser.parse_args()

    # parse stdout and read data into numpy arrays
    brick= args.logfile.split('.')[-1]
    print('brick=',brick)
    fmem= os.path.join(args.outdir,'mem_%s_rs%s.txt' % (brick,args.rsdir))
    ftime= os.path.join(args.outdir,'time_%s_rs%s.txt' % (brick,args.rsdir))
    # multi node
    # grep "runbrick.py starting at" bb_multi.o2907746
    # grep "Stage writecat finished:" bb_multi.o2907746
    if os.path.exists(fmem) and os.path.exists(ftime):
        print('using existing files:\n%s\n%s' % (fmem,ftime))
    else:
        bash_result("grep 'Resources for' %s -A 2|grep -e 'Wall:' -e 'Resources' > %s" % \
                    (args.logfile,fmem))
        bash_result("grep -e 'Resources for stage' -e 'Total serial Wall' -e 'Total parallel Wall' -e 'Grand total Wall' -e 'Grand total CPU utilization' %s > %s" % \
                    (args.logfile,ftime))
    mem=parse_tractor_profile(fmem)
    tm=parse_tractor_profile(ftime)
    raise ValueError
    # plots
    ncores= str(bash_result("grep 'Command-line args:' %s|cut -d ',' -f 11|tail -n 1" % args.logfile) )
    ncores= int( ncores.replace('"','').replace("'",'') )
    # write timing info to text file and create a pickle file
    fout=open('result_%s_rs%s.txt' % (brick,args.rsdir),'a')
    fout.write('#nodes cores times[sec]:tims fitblobs total\n')
    fout.write('1 %d %.2f %.2f %.2f\n' % (ncores,tm['total'][0],tm['total'][4],tm['total'].sum()))
    fout.close()

##fn="$1"
#grep -e "Resources for" ${fn} -A 2|grep -e "Wall:" -e "Resources" > ${fn}_mem.txt
#grep -e "Resources for stage" -e "Total serial Wall" -e "Total parallel Wall" -e "Grand total Wall" -e "Grand total CPU utilization" ${fn} > ${fn}_time.txt
#echo 'done!'
