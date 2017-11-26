"""
Generally useful functions for other modules or repos
"""
import matplotlib
import matplotlib.pyplot as plt
import os
import pandas as pd

def inJupyter():
    return 'inline' in matplotlib.get_backend()
    
def save_png(outdir,fig_id, tight=True):
    path= os.path.join(outdir,fig_id + ".png")
    if not os.path.isdir(outdir):
        os.makedirs(dirname)
    print("Saving figure", path)
    if tight:
        plt.tight_layout()
    plt.savefig(path, format='png', dpi=150)
    #plt.savefig(path, format='png',box_extra_artists=[xlab,ylab],
    #            bbox_inches='tight',dpi=150)
    if not inJupyter():
        plt.close()

def dobash(cmd):
    print('UNIX cmd: %s' % cmd)
    if os.system(cmd): raise ValueError

def fits2pandas(tab,attrs=None):
    """converts a fits_table into a pandas DataFrame

    Args:
        tab: fits_table()
        attrs: attributes or column names want in the DF
    """
    d={}
    if attrs is None:
        attrs= tab.get_columns()
    for col in attrs:
        d[col]= tab.get(col)
    df= pd.DataFrame(d)
    # Fix byte ordering from fits
    # https://stackoverflow.com/questions/18599579/pulling-multiple-non-consecutive-index-values-from-a-pandas-dataframe
    df= df.apply(lambda x: x.values.byteswap().newbyteorder())
    return df

def get_brickdir(outdir,obj,brick):
    return os.path.join(outdir,obj,brick[:3],brick)

def get_outdir_runbrick(outdir,obj,brick,rowstart,
                do_skipids='no',do_more='yes'):
    """diretory obiwan/runbrick will write results to

    Returns path to like outdir/obj/bri/brick/rs0
    """
    # Either rs or skip_rs
    if do_skipids == 'no':
        final_dir= "rs%s" % str(rowstart)
    elif do_skipids == 'yes':
        final_dir= "skip_rs%s" % str(rowstart) 
    # if specified minimum id, running more randoms
    if do_more == 'yes':
        final_dir= "more_"+final_dir
    #return os.path.join(get_brickdir(outdir,obj,brick),
    #                  final_dir)
    return os.path.join(outdir,obj,brick[:3],brick,final_dir)

def get_outdir_final(outdir,obj,brick,rowstart,
                do_skipids='no',do_more='yes'):
    """directory obiwan/runbrick outputs will be repackaged into

    Returns path to like outdir/obj/coadd/bri/brick/rs0
    """
    # Either rs or skip_rs
    if do_skipids == 'no':
        final_dir= "rs%s" % str(rowstart)
    elif do_skipids == 'yes':
        final_dir= "skip_rs%s" % str(rowstart) 
    # if specified minimum id, running more randoms
    if do_more == 'yes':
        final_dir= "more_"+final_dir
    #return os.path.join(get_brickdir(outdir,obj,brick),
    #                  final_dir)
    return os.path.join(outdir,obj,brick[:3],brick,final_dir)
