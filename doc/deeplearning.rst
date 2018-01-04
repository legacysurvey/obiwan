.. raw:: html
  
  <h1>Deep Learning</h1>

Can a CNN tell the difference between fake galaxies (that I inject into real imaging data) and real ones?

.. contents:: Table of Contents
  :depth: 2


Background
------------

I am getting my PhD in astrophysics from UC Berkeley. My research focuses on observational cosmology, which for me boils down to analyzing images of the night sky. Specifically, building an pipeline that identifies and fits models to stars and galaxies in multi-wavelength images. From the location and properties of the galaxies it is possible to answer questions like, “how has the expansion rate of the Universe changed over time?”. The `Sloan Digital Sky Survey <http://www.sdss.org>`_ (SDSS) showed that you simply need a large enough sample of galaxies over a large enough region of the sky. They imaged 1/3 of the night sky to detect about 1 million galaxues up to 6.7 billion light years away. My team is running the `Legacy Survey <http://www.legacysurvey.org>`_, which will detect 30x more galaxies than SDSS and up to 6x fainter ones (an additional 3.2 billion lightyears away). 

Our goal is a 2D map showing the positions of all the galaxies from about 100 TBs of images. This is scientifically important because we can measure the distance to each galaxy from it’s spectrum (e.g. how bright it is at each wavelength of visible light). Distance makes our map 3D and from this we can measure the expansion rate of the Universe at different points in time (`Eisenstein & Hu 1998 <https://arxiv.org/abs/astro-ph/9709112>`_; `Eisenstein et al. 2005 <https://arxiv.org/abs/astro-ph/0501171>`_; `Seo & Eisenstein 2007 <https://arxiv.org/abs/astro-ph/0701079>`_; `Butler et al. 2017  <https://arxiv.org/abs/1607.03150>`_).

The primary goal of my thesis is to measure the statistical bias and variance of the `pipeline <https://github.com/legacysurvey/legacypipe>`_ responsible for detecting and modeling galaxies and stars in the images. How does our completeness depend on whether a galaxy is bright or faint, blue or red, big or small? How well does our pipeline handle image artifacts, instrument issues, or transient objects?

To answer all of this, John Moustakas, a professor at Siena College, and I created the `obiwan code <https://github.com/legacysurvey/obiwan>`_ to do Monte Carlo simulations of our `pipeline <https://github.com/legacysurvey/legacypipe>`_. I inject fake (but realistic looking) galaxies into the current data set of images, at random locations, and then rerun the pipeline on the modified images.

The success of all this depends on whether or not the fake galaxies are representative of the true galaxy population. If not, I will measure our pipeline’s bias and variance for a sample of objects that do not exist in the data. I need to show that the fake galaxies have the same properties as  the real ones and that they "look" the same. Enter the Convolutional Neural Network (CNN)...

The Problem
------------

This is a supervised binary classification problem of labeled multi-wavelength images. The training set is at least 60 million examples. Each example is six images (64 x 64 x 6 pixels): one image for each of the three wavelength bands and three images to encode camera artifacts and the variance of each pixel. 

I am training a CNN to predict whether an image of a galaxy is model-generated (fake).  

Novel 
------

The project is novel because I want the CNN to do *poorly* even when well trained. If the CNN cannot do better than guessing the most numerous class, then the fake galaxies "look like" real ones and my model is representative of reality.

This has attracted the interest of other teammates because I am testing how realistic a model is, by injecting model-generated images into real data, and comparing against real examples elsewhere in the data. This generalizes to any model for a galaxy, star, or image artifact the group may have. 


Examples of Training Data
--------------------------

The galaxies below are some of the brightest in the training set (99th percentile in brightness). :ref:`More Examples <more-examples>` of galaxies are at the end. Those are fainter having 75th, 50th, 25th, and 1st percentiles in brightness.

.. figure:: _static/fake_real_mosaic_istart_0.png
   :width: 600 px
   :figwidth: 600 px
   :align: center

   **Figure 1.** The label for each image is on the left (R for Real and F for Fake)  and its corresponding g-band magnitude is the number on the right (the smaller the number, the brighter the galaxy). Each row represents a single galaxy imaged at three different wavelengths. The color-image (left most panel) shows the colors you would see by eye, while the black and white-images (right six panels) are the training data of individual wavelength (g, r, z) and corresponding inverse variance (ivar g, r, z) images. Finally, consecutive rows of R and F (rows 1 and 2, 3 and 4, etc.) have similar g-band magnitudes so that a fair comparison can be made.

These examples and those :ref:`at the end <more-examples>` reveal at least two challenges for the CNN.

#. Only the central-object matters, but there are many off-center objects in the images. These are random background sources, often bright galaxies or stars that we are not interested in. 
#. These galaxies are very faint. The CNN must be able to dig out the low Signal to Noise sources.

CNN Architecture
-----------------

As a starting point, I used TensorFlow to build a CNN similar to LeNet-5 with the following architecture: 

.. list-table:: 
   :widths: auto
   :header-rows: 1
   :align: left

   * - Layer
     - Feature Maps
     - Size
     - Kernel Size
     - Stride
     - Activation Function
   * - Input Image
     - 6
     - 64x64
     - 
     - 
     - 
   * - Convolution
     - 18
     - 64x64
     - 7x7
     - 1
     - ReLU
   * - Avg. Pooling
     - 18
     - 32x32
     - 7x7
     - 2
     - ReLU
   * - Convolution
     - 36
     - 32x32
     - 7x7
     - 1
     - ReLU
   * - Avg. Pooling
     - 36
     - 16x16
     - 7x7
     - 2
     - ReLU
   * - Convolution
     - 54
     - 16x16
     - 7x7
     - 1
     - ReLU
   * - Avg. Pooling
     - 54
     - 8x8
     - 7x7
     - 2
     - ReLU
   * - Fully Connected
     - 
     - 64 
     - 
     - 
     - ReLU
   * - Fully Connected
     - 
     - 2 
     - 
     - 
     - Softmax

The input image has 64 x 64 x 6 pixels. With three convolution/pooling layers, the CNN is much shallower than the ImageNet ILSVRC winners, so in addition to tuning the number of feature maps, kernel size, stride, etc., I plan to make it deeper.

Intel Xeon Phi
-----------------------------------

I have used the Cray XC30 (`Edison <http://www.nersc.gov/users/computational-systems/edison/>`_) and Cray XC40 (`Cori <http://www.nersc.gov/users/computational-systems/cori/>`_) supercomputers at the National Energy Research Scientific Computing Center (NERSC) for the majority of my thesis work. With almost 10,000 Intel Xeon Phi processor nodes on Cori, NERSC Staff are particularly interested in helping users optimize their codes for Xeon Phi. 

I decided to train on Xeon Phi when NERSC/Intel released optimized installs of many of the popular machine learning libraries (Caffe, TensorFlow, Theano, Torch, see `full list <http://www.nersc.gov/users/data-analytics/data-analytics-2/deep-learning/using-tensorflow-at-nersc>`_). I created an initial training set of 2048 images with an equal number of fake and real galaxies. The images are float32 so I stored every 512 examples in a file, thinking that a 50 MB file would fit in memory of most machines.

It takes about 3 minutes to train 4 epochs of 2048 images on a single Xeon Phi node (68 hardware cores). For hundreds of nodes, I plan on training on a different minibatch with each MPI task, updating a global set of weights once all MPI tasks finish, then repeating. Although NERSC Staff are scaling ResNet-50 and DCGAN to 1024 Xeon Phi nodes, multi-node support is not yet available to users.

Fortunately, the NERSC Staff have volunteered my CNN for non-benchmark multi-node testing. I hope to begin multi-node training soon. 


TensorBoard & Profiling
------------------------

The accuracy, loss, and graph for the 4 epochs of training on 2048 images is shown with TensorBoard, below. The different colors correspond to me restarting the training twice to demonstrate that the checkpoints are working.

.. figure:: _static/tensorboard_scalars.png
   :width: 75 %
   :figwidth: 75 %
   :align: center

   Accuracy and loss with TensorBoard 

.. figure:: _static/tensorboard_graph.png
   :width: 75 %
   :figwidth: 75 %
   :align: center

   Graph with TensorBoard

I also profile my CNN using TensorFlow's `timeline <https://stackoverflow.com/questions/34293714/can-i-measure-the-execution-time-of-individual-operations-with-tensorflow>`_ object. This times each node of the graph and writes a json file. Google Chrome automatically displays the profiling info by going to `chrome://tracing`, clicking `load`, and selecting the file. It looks like this for the 4 training epochs:

.. figure:: _static/prof_tensorflow.png
   :width: 90 %
   :figwidth: 90 %
   :align: center

   Profiling with TensorFlow's `timeline` and Google Chrome 


.. _deep-learn-instructions:

Instructions
-----------------

These are the instructions are for creating the training set and training the CNN at NERSC.

Create Training Data
"""""""""""""""""""""""""
* :mod:`obiwan.dplearn.create_training` (source code: `create_training.py <https://github.com/legacysurvey/obiwan/blob/master/py/obiwan/dplearn/create_training.py>`_) saves 64x64 pixels cutouts of each source in a Data Release to an HDF5 file, indexed by its unique tractor id. One HDF5 file per brick. This is done for real galaxies using an official Data Release and again for fake ones using the results from Obiwan. 

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
"""""""""""""""""""""""

* :mod:`obiwan.dplearn.split_testtrain` (source code: `split_testtrain.py <https://github.com/legacysurvey/obiwan/blob/master/py/obiwan/dplearn/split_testtrain.py>`_) randomly shuffles the real and fake galaxies in the above HDF5 files, does a 80% training/20% test split, and repackages the results in numpy binary files.

It uses mpi4py so the same SLURM job can be used, expect with::

    srun -n ${tasks} -c 1 python split_testtrain.py \
         --bricks_fn bricks.txt --nproc ${tasks} \
         --real_dir /global/cscratch1/sd/kaylanb/obiwan_out/dr5_hdf5 \
         --sim_dir /global/cscratch1/sd/kaylanb/obiwan_out/elg_dr5_coadds \
         --save_dir /global/cscratch1/sd/kaylanb/obiwan_out/dr5_testtrain

The resulting numpy files are on at NERSC:
* /global/cscratch1/sd/kaylanb/obiwan_out/dr5_testtrain

The training data are named `[xy]train_[0-9]+.npy` and have 512 `64x64x6` examples per file. The test data are named `[xy]test_[0-9]+.npy`.


Train the CNN
"""""""""""""""""

* :mod:`obiwan.dplearn.cnn` (source code: `cnn.py <https://github.com/legacysurvey/obiwan/blob/master/py/obiwan/dplearn/cnn.py>`_) trains the CNN using TensorFlow. The following SLURM job will run on a single Knights Landing (KNL) node using 68 threads ("srun" is not needed because this is a single node job)::

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

.. _more-examples:

More Examples
--------------------------------

Galaxies with *75th* perentile in brightness:

.. figure:: _static/fake_real_mosaic_istart_64.png
   :width: 600 px
   :figwidth: 600 px
   :align: center

Galaxies with *50th* perentile in brightness:

.. figure:: _static/fake_real_mosaic_istart_112.png
   :width: 600 px
   :figwidth: 600 px
   :align: center

Galaxies with *25th* perentile in brightness:

.. figure:: _static/fake_real_mosaic_istart_208.png
   :width: 600 px
   :figwidth: 600 px
   :align: center

Galaxies with *1st* percentile in brightness (some of the *faintest* galaxies in the training set):

.. figure:: _static/fake_real_mosaic_istart_254.png
   :width: 600 px
   :figwidth: 600 px
   :align: center

