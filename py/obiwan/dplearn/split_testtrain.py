"""
Shuffles the fake and real galaxies in the HDF5, does train/test split, repackages as numpy binaries
"""

import numpy as np
import os
from glob import glob
import h5py

def get_data(f,keys):
    """Returns numpy array of shape [len(keys),64,64,6]
    
    Args:
        f: hdf5 file object
    """
    return np.array([np.stack([f[key+'/img'],f[key+'/ivar']],axis=-1).reshape((64,64,6))
                     for key in keys])

def write_traintest(brick,real_dir,sim_dir,save_dir,
                    n_train=256,n_test=64):
    """Writes xtrain_1.npy,xtest_1.npy,... for a given brick

    Args:
        brick: brickname
        real_dir: path to hdr5 dir for real galaxies
        sim_dir: ... simulated galaxies
        save_dir: where to write the bri/brick/xtrain.npy, ..., files
    """
    bri=brick[:3]
    f_real= h5py.File(os.path.join(real_dir,'hdf5',bri,brick,
                                   'img_ivar_grz.hdf5'),
                      'r')
    f_sim= h5py.File(os.path.join(sim_dir,'hdf5',bri,brick,
                                   'img_ivar_grz.hdf5'),
                      'r')
    keys_real= np.array(list(f_real.keys()))
    keys_sim= np.array(list(f_sim.keys()))
    n_samples= np.min([keys_real.size,keys_sim.size])
    print('n_samples=',n_samples)

    # Shuffle real
    ind=np.arange(keys_real.size).astype(int)
    np.random.shuffle(ind)
    keys_real= keys_real[ind]

    # Sort sim by id is equiv to shuffling
    keys_sim= np.sort(keys_sim.astype(np.int32)).astype(str)

    # Take equal size samples
    keys_real= keys_real[:n_samples]
    keys_sim= keys_sim[:n_samples]

    chunk_size= n_train + n_test
    print('chunk_size=',chunk_size)
    # chunks of equal size exept last chunk is whateve is left over
    #ind= np.array_split(np.arange(n_samples),n_samples // chunk_size + 1)
    #print('ind.shape=',ind[0].shape)
    for cnt in range(n_samples // chunk_size + 1): # in enumerate(ind):
        slc= slice(cnt*chunk_size,(cnt+1)*chunk_size)
        print('slice= ',slc)
        print('Brick=%s, chunk=%d' % (brick,cnt+1)) 
        # Shape [chunk_size,64,64,6]
        Xreal= get_data(f_real,keys_real[slc])
        Yreal= np.zeros(len(keys_real[slc]))
        Xsim= get_data(f_sim,keys_sim[slc])
        Ysim= np.ones(len(keys_sim[slc]))
        print(Xreal.shape,Xsim.shape,Yreal.shape,Ysim.shape)
        # Combine and shuffle
        x= np.vstack([Xreal,Xsim])
        y= np.hstack([Yreal,Ysim])
        print(x.shape,y.shape)
        shuff= np.arange(x.shape[0])
        np.random.shuffle(shuff)
        x= x[shuff,...]
        y= y[shuff]
        # Split
        isplit= 2*n_train
        if x.shape[0] < 2*chunk_size:
            isplit= int(float(n_train/chunk_size) * x.shape[0])
        print('isplit=',isplit)
        Xtrain,Xtest= x[:isplit,...],x[isplit:,...]
        Ytrain,Ytest= y[:isplit],y[isplit:]
        print(Xtrain.shape,Xtest.shape,Ytrain.shape,Ytest.shape)
        # Save
        dr= os.path.join(save_dir,'testtrain',bri,brick)
        try: 
            os.makedirs(dr)
        except OSError:
            print('directory already exists: ',dr)
        for data,name in zip([Xtrain,Xtest,Ytrain,Ytest],
                             ['xtrain','xtest','ytrain','ytest']):
            fn=os.path.join(dr,'%s_%d.npy' % (name,cnt+1))
            np.save(fn,data)
            print('Wrote %s' % fn)



if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--nproc', type=int, default=1, help='set to > 1 to run mpi4py') 
    parser.add_argument('--bricks_fn', type=str, default=None, help='specify a fn listing bricks to run, or a single default brick will be ran') 
    parser.add_argument('--real_dir', type=str, required=True, help='path to hdr5 dir for real galaxies') 
    parser.add_argument('--sim_dir', type=str, required=True, help='ditto for sim galaxies') 
    parser.add_argument('--save_dir', type=str, required=True, help='where to write the bri/brick/xtrain.npy, ..., files') 
    args = parser.parse_args()

    if args.bricks_fn:
        bricks= np.loadtxt(args.bricks_fn,dtype=str)
    else:
        bricks=['1211p060']

    if args.nproc > 1:
        from mpi4py.MPI import COMM_WORLD as comm
        bricks= np.array_split(bricks, comm.size)[comm.rank]

    for brick in bricks:
        write_traintest(brick,args.real_dir,args.sim_dir,args.save_dir)


