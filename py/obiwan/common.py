"""
Generally useful functions for other modules or repos
"""
import matplotlib
import matplotlib.pyplot as plt
import os

def inJupyter():
    return 'inline' in matplotlib.get_backend()
    
def save_png(outdir,fig_id):
    path= os.path.join(outdir,fig_id + ".png")
    if not os.path.isdir(outdir):
        os.makedirs(dirname)
    print("Saving figure", path)
    plt.tight_layout()
    plt.savefig(path, format='png', dpi=150)
    #plt.savefig(path, format='png',box_extra_artists=[xlab,ylab],
    #            bbox_inches='tight',dpi=150)
    if not inJupyter():
        plt.close()

def dobash(cmd):
  print('UNIX cmd: %s' % cmd)
  if os.system(cmd): raise ValueError

def get_brickdir(outdir,obj,brick):
  return os.path.join(outdir,obj,brick[:3],brick)

def get_savedir(outdir,obj,brick,rowstart,
                do_skipids='no',do_more='yes'):
  # Either rs or skip_rs
  if do_skipids == 'no':
    final_dir= "rs%s" % str(rowstart)
  elif do_skipids == 'yes':
    final_dir= "skip_rs%s" % str(rowstart) 
  # if specified minimum id, running more randoms
  if do_more == 'yes':
    final_dir= "more_"+final_dir
  return os.path.join(get_brickdir(outdir,obj,brick),
                      final_dir)


