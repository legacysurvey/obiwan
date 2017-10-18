if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg') #display backend
from argparse import ArgumentParser
import numpy as np
import matplotlib.pyplot as plt
import fitsio
import glob
import os
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from astropy import units
from astropy.coordinates import SkyCoord

from astrometry.util.fits import fits_table, merge_tables
from tractor.brightness import NanoMaggies
from theValidator.catalogues import ConvertTractor,CatalogueFuncs, Matcher
from theValidator.plots import Kaylans,Dustins
#import plots_OldDataModel as plots

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('--ref_cat_list',action='store',required=True,help="text file listing refernce tractor catalogues, one per line")
    parser.add_argument('--obs_cat_list',action='store',required=True,help="text file listing tractor catalogues to compare with, one per line")
    parser.add_argument('--ref_name',action='store',default='REF',required=False,help="name for reference catalogues")
    parser.add_argument('--obs_name',action='store',default='OBS',required=False,help="name for catalogues comparing to")
    parser.add_argument('--shuffle',type=int,default=None,required=False,help="set to an integer to randomly read of the N catalogues instead of all N")
    args = parser.parse_args()

    ref= CatalogueFuncs().stack(args.ref_cat_list, shuffle=args.shuffle)
    obs= CatalogueFuncs().stack(args.obs_cat_list, shuffle=args.shuffle)

    # New naming scheme
    conv= ConvertTractor()
    if conv.isOld(ref.get_columns()):
        print('Converting ref')
        ref= conv.old_to_new(ref)
    if conv.isOld(obs.get_columns()):
        print('Converting obs')
        obs= conv.old_to_new(obs)    
    
    print('matching')
    imatch,imiss,d2d= Matcher().match_within(ref,obs) #,dist=1./3600)

    print('adding columns')
    #CatalogueFuncs().set_mags_OldDataModel(ref)
    CatalogueFuncs().set_mags(ref)
    CatalogueFuncs().set_mags(obs)
            
    # Plots 
    oldformat_any=False
    #if hasattr(ref,'decam_flux') or hasattr(obs,'decam_flux'):
    #    oldformat_any=True 
    k=Kaylans(ref,obs,imatch,
              ref_name=args.ref_name,obs_name=args.obs_name,savefig=True,
              oldformat_any=oldformat_any)
    d=Dustins(ref,obs,imatch,d2d,
              ref_name=args.ref_name,obs_name=args.obs_name,savefig=True,
              oldformat_any=oldformat_any)
