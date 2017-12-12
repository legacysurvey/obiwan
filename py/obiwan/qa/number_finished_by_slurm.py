import pandas as pd
import os
import re

def add_fits(text):
    return text+'.fits'

def trac_fns(slurm_fn):
    with open(slurm_fn,'r') as foo:
        text= foo.read()
    return (pd.Series(re.findall(r'Logging to:.*?\n',text))
             .str.replace(r'Logging to:\s','')
             .str.strip()
             .str.replace('logs','tractor')
             .str.replace(r'log\.','tractor-')
             .apply(add_fits)
            ).values

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser= ArgumentParser()
    parser.add_argument('--slurm_fn', type=str, required=True)
    args = parser.parse_args()

    tractor_fns= trac_fns(args.slurm_fn)
    cnts= [1 if os.path.exists(fn) else 0 
           for fn in tractor_fns]
    print('1: tractor.fits exists, 0: doesnt')
    print(pd.Series(cnts).value_counts())
