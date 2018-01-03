.. raw:: html
  
  <h1>Training Instructions</h1>

The following instructions are for training the CNN described on the :doc:`deep learning page <../deeplearning>`.

.. contents:: Table of Contents
  :depth: 2


Create Training Data
^^^^^^^^^^^^^^^^^^^^^^^^^

`create_training.py <https://github.com/legacysurvey/obiwan/blob/master/py/obiwan/dplearn/create_training.py>`_ creates 64x64 pixels cutouts for each source in the coadd/ directory of a Data Release, and saves the cutouts (indexed by unique tractor id) in an HDF5 file. One HDF5 file is written per brick. This is done for real galaxies using an official Data Release and again for fake ones using the results from Obiwan. 

Fake galaxies occupy the narrow region of parameter space we are interested in, while real galaxies do not. The only difference in procedure between building the fake and real training sets is removing real galaxies that are outside the parameter space of interest. 

There are millions of fake and real galaxy images, so the script uses mpi4py and scales well to a few hundred Haswell nodes. I created about 1 million **real** galaxy examples using 50 Haswell nodes for 1 hour with the following SLURM job script::

    #SBATCH -p regular
    #SBATCH -N 50
    #SBATCH -t 01:00:00
    #SBATCH --account=desi
    #SBATCH -J train
    #SBATCH -L SCRATCH,project
    #SBATCH -C haswell

    let tasks=32*${SLURM_JOB_NUM_NODES}

    # NERSC / Cray / Cori / Cori KNL things
    export KMP_AFFINITY=disabled
    export MPICH_GNI_FORK_MODE=FULLCOPY
    export MKL_NUM_THREADS=1
    export OMP_NUM_THREADS=1

    srun -n ${tasks} -c 1 python create_training.py \
         --which tractor --bricks_fn bricks.txt --nproc ${tasks}
         --savedir /global/cscratch1/sd/kaylanb/obiwan_out/dr5_hdf5

For **fake** galaxies, simply replace "--which tractor" with "--which sim". The resulting HDF5 files are on at NERSC:

* real from DR5: /global/cscratch1/sd/kaylanb/obiwan_out/dr5_hdf5
* fake from Obiwan using DR5: /global/cscratch1/sd/kaylanb/obiwan_out/elg_dr5_coadds/hdf5 


Split Train/Test
^^^^^^^^^^^^^^^^^

`split_testtrain.py <https://github.com/legacysurvey/obiwan/blob/master/py/obiwan/dplearn/split_testtrain.py>`_ randomly shuffles the real and fake galaxies in the above HDF5 files, does a 80% training/20% test split, and repackages the results in numpy binary files.

Again using mpi4py, the same SLURM job can be used. This time with::

    srun -n ${tasks} -c 1 python split_testtrain.py \
         --bricks_fn bricks.txt --nproc ${tasks} \
         --real_dir /global/cscratch1/sd/kaylanb/obiwan_out/dr5_hdf5 \
         --sim_dir /global/cscratch1/sd/kaylanb/obiwan_out/elg_dr5_coadds/hdf5 \
         --save_dir /global/cscratch1/sd/kaylanb/obiwan_out/dr5_testtrain

The resulting numpy files are on at NERSC:
* /global/cscratch1/sd/kaylanb/obiwan_out/dr5_testtrain

The training data are named `[xy]train_[0-9]+.npy` and have 512 `64x64x6` examples per file. The test data are named `[xy]test_[0-9]+.npy`.


Train the CNN
^^^^^^^^^^^^^^^

`cnn.py <https://github.com/legacysurvey/obiwan/blob/master/py/obiwan/dplearn/cnn.py>`_ trains the CNN using TensorFlow. The following SLURM job will run on a single Knights Landing (KNL) node using 68 threads ("srun" is not needed because this is a single node job)::

    #!/bin/bash
    #SBATCH -N 1
    #SBATCH -C knl,quad,cache
    #SBATCH -p debug
    #SBATCH -J tf
    #SBATCH -t 00:30:00

    module load tensorflow/intel-head
    export OMP_NUM_THREADS=68
    export KMP_AFFINITY="granularity=fine,verbose,compact,1,0"
    export KMP_SETTINGS=1
    export KMP_BLOCKTIME=1
    export isKNL=yes

    python cnn.py --outdir /global/cscratch1/sd/kaylanb/obiwan_out/cnn

This will write three sets of metadata:

* checkpoints: /global/cscratch1/sd/kaylanb/obiwan_out/cnn/**ckpts**
* tensorboard logs: /global/cscratch1/sd/kaylanb/obiwan_out/cnn/**logs**
* profiling info: /global/cscratch1/sd/kaylanb/obiwan_out/cnn/**prof**

If checkpoints files exists, the CNN will restart training from there and the appropriate epoch and batch will be selected.

TensorBoard
^^^^^^^^^^^^^

You can use `ssh` to launch TensorBoard on your laptop but using the log files on NERSC. Following NERSC's `instructions <http://www.nersc.gov/users/data-analytics/data-analytics-2/deep-learning/using-tensorflow-at-nersc/#TensorBoard>`_, do::

    ssh <user>@cori.nersc.gov
    module load tensorflow/intel-head
    tensorboard --logdir=/global/cscratch1/sd/kaylanb/obiwan_out/cnn/logs --port 9998
    # new terminal
    ssh -L 9998:localhost:9998 <user>@cori.nersc.gov

Then paste `http://0.0.0.0:9998/` in your browser. It should look like this

.. figure:: ../_static/fake_real_mosaic_istart_0.png
   :width: 100 %
   :figwidth: 100 %
   :align: center

   Tensorboard example

Profile
^^^^^^^^

We profile our CNN with TensorFlow's `timeline <https://stackoverflow.com/questions/34293714/can-i-measure-the-execution-time-of-individual-operations-with-tensorflow>`_ object, which times each node in the graph. The times are saved to a json file "timing.json"

You can inspect the file with your Chrome browser::

    scp <user>@cori.nersc.gov:/global/cscratch1/sd/kaylanb/obiwan_out/cnn/prof/timing.json ./

Then go to url `chrome://tracing`, click `load`, and select "timing.json". It should look like this

.. figure:: ../_static/prof_tensorflow.png
   :width: 100 %
   :figwidth: 100 %
   :align: center

   Profiling example 





