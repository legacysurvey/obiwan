.. raw:: html
  
  <h1>Deep Learning</h1>

Can a CNN tell the difference between fake galaxies (that I inject into real imaging data) and real ones?

.. contents:: Table of Contents
  :depth: 2


Background
^^^^^^^^^^^^

I am getting my PhD in astrophysics from UC Berkeley. My research focuses on observational cosmology, which for me, means building an image analysis pipeline that identifies and measures properties of stars and galaxies in multi-wavelength images of the night sky. It is possible (believe it or not!) to answer questions like, “how has the expansion rate of the Universe changed over time?”, using a statistical survey of galaxies over a large enough volume of the Universe. The `Sloan Digital Sky Survey <http://www.sdss.org>`_ (SDSS) carried out such a survey by imaging about one-third of the night sky and detecting galaxies up to 6.7 billion light years away. My team is leading the `Legacy Survey <http://www.legacysurvey.org>`_. It is similar to SDSS but we use three telescopes instead of one and will detect about 30x more galaxies that can be up to 6x fainter and 3.2 billion light years further away than the SDSS galaxies. Our goal is a 2D map showing the positions of 30 million galaxies from about 100 TBs of imaging data. This is scientifically important because we can measure the distance to each galaxy from it’s spectrum (e.g. how bright it is at each wavelength of visible light). From our 3D map (distance is the third dimension), we can measure the expansion rate of the Universe at different points in time (`Eisenstein & Hu 1998 <https://arxiv.org/abs/astro-ph/9709112>`_; `Eisenstein et al. 2005 <https://arxiv.org/abs/astro-ph/0501171>`_; `Seo & Eisenstein 2007 <https://arxiv.org/abs/astro-ph/0701079>`_; `Butler et al. 2017  <https://arxiv.org/abs/1607.03150>`_).

Our survey is one of the first cosmological surveys to be entirely open source. Anyone can `download <http://archive.noao.edu/search/query>`_ the raw and calibrated images from all three telescopes within a few days of being observed. Every six months, we publicly release the results of our open source `pipeline <https://github.com/legacysurvey/legacypipe>`_ that detects and models all of the galaxies and stars in the imaging data set. I have developed many parts of the pipeline, and successfully ran it on about 30 TBs of imaging to produce our `4th public release <http://legacysurvey.org/dr4/description>`_ in June 2017. 

The primary goal of my thesis is to measure the statistical bias and variance of our pipeline. For instance, how does our survey completeness depend on whether a galaxy is bright or faint, blue or red, big or small? How well does our pipeline handle the various artifacts, every evolving list of instrument issues, and transient sources that are in our images?

To answer all of this, a professor at Siena College, John Moustakas, and I created the `obiwan code <https://github.com/legacysurvey/obiwan>`_ to do Monte Carlo simulations of our `pipeline <https://github.com/legacysurvey/legacypipe>`_. We inject fake (but realistic looking) galaxies into the current data set of images, at random locations, and then rerun the pipeline on the modified images. We expect to model about 30 million galaxies in the final data set, so we need about 10x as many fake objects to reduce standard errors on the measurements. That’s about 300 million fake galaxies.

The success of all this depends on whether or not the fake galaxies are representative of the true galaxy population. If not, we will have spent 5-10 million CPU hours measuring our pipeline’s bias and variance for a sample of objects that does not generalize to our actual sample of real galaxies. The fake and real galaxies have similar age, brightness, color, and shape distributions; however, these metrics do not necessarily say whether the fake and real galaxies “look” the same. 

There are 30 million real galaxies and at least that many fake ones. Each galaxy is really six images (not one): one image for each of the three wavelength bands and an image for each of those that encodes camera artifacts and the variance of each pixel. I need an algorithm that can visually “look at” all of the images, so a Convolutional Neural Network (CNN) seems like the obvious choice.


Novel 
^^^^^^

This Deep Learning project is novel because of the training set. I am testing how realistic my model for a specific type of galaxy is, by injecting model-generated images of that galaxy into the imaging data and comparing against real examples of that galaxy. I can create an unlimited number of model-generated “fake” galaxies, so the size of the training set is limited by the number of “real” galaxies, about 30 million.

This project has also attracted the interest on other teammates because it is a great test-bed for user-defined models. Given a model for any type of object (galaxy, star, or image artifact), it can be validated against real examples of that object.


CNN Architecture
^^^^^^^^^^^^^^^^^

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

Source Code
^^^^^^^^^^^^

Follow these :doc:`instructions <howto/deeplearn>` for training the CNN.

* :mod:`obiwan.dplearn.create_training` creates the training data. Source code: `create_training.py <https://github.com/legacysurvey/obiwan/blob/master/py/obiwan/dplearn/create_training.py>`_. 
* :mod:`obiwan.dplearn.cnn` trains the CNN. Source code: `cnn.py <https://github.com/legacysurvey/obiwan/blob/master/py/obiwan/dplearn/cnn.py>`_.

Training Data
^^^^^^^^^^^^^

This is a supervised binary classification problem using labeled multi-wavelength imaging data. There are at least 30 million examples for each of the two classes, e.g. “real” or “fake” galaxy, and each example has 64 x 64 x 6 pixels. “64 x 64” because this is the smallest size that encloses all the galaxy examples and “x 6” because of the multi wavelength data. My question is, “do fake galaxies (which I create and inject into the imaging data) look like real galaxies (the ones actually in the data)?” If so, the CNN cannot do better than guessing the most numerous class.

The *brightest* galaxies (99th percentile in brightness) are probably the easiest to classify, so let's look at a subset of these first: 

.. figure:: _static/fake_real_mosaic_istart_0.png
   :width: 600 px
   :figwidth: 600 px
   :align: center

   **Figure 1.** The label for each image is on the left (R for Real and F for Fake)  and its corresponding g-band magnitude is the number on the right (the smaller the number, the brighter the galaxy). Each row represents a single galaxy imaged at three different wavelengths. The color-image (left most panel) shows the colors you would see by eye, while the black and white-images (right six panels) are the training data of individual wavelength (g, r, z) and corresponding inverse variance (ivar g, r, z) images.

The :ref:`Additional Training Examples <additional-examples>` include four more panels of fainter galaxies: the 75th, 50th, 25th, 1st percentiles in brightness.

To compare similar R and F images, consecutive rows have similar g-band magnitudes. For example, rows 1 and 2 have nearly the same g-band magnitude, and the same for rows 3 and 4, etc.

Challenges
^^^^^^^^^^

There are at least two challenges apparent in the training examples above. First, there are many off-center objects in the images. These are random background sources, often bright galaxies or stars that we are not interested in. The CNN must learn to just care about the central-object. Second, the galaxies in our sample are very faint, so even without background sources the CNN must be able to detect a low Signal to Noise (S/N) central-galaxy.


Training the CNN
^^^^^^^^^^^^^^^^^^

I created an initial data set of two million images with an equal number of “fake” and “real” examples. I randomly split that into 80% training and 20% test, storing every 512 examples (32 bit floating point) in a numpy binary file so this 50 MB file would fit in memory on most machines.

I am currently training my CNN with a batch size of 16 on Xeon Phi (Knights Landing, KNL) CPUs. This non-GPU choice was motivated by the recent addition of KNL nodes to the National Energy Research Scientific Computing Center’s (NERSC) Cray XC40 supercomputer “Cori”, and the opportunity for NERSC users to see how well their codes can scale on the new system. 

NERSC has installed many of the popular machine learning packages (Caffe, TensorFlow, Theano, Torch, see `full list <http://www.nersc.gov/users/data-analytics/data-analytics-2/deep-learning/using-tensorflow-at-nersc>`_) on Cori and optimized them for KNL. I can only train on 1 node (68 threads) because multi-node support is “coming soon,” but I’ve been told that I should be able to begin multi-node training soon because they can now scale ResNet-50 and DCGAN to 1024 KNL nodes. When that happens, I plan to assign a different batch to each MPI task, update a global set of weights after each back propagation step, and repeat.

.. _additional-examples:

Additional Training Examples
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A subset of galaxies with *75th perentile* in brightness:

.. figure:: _static/fake_real_mosaic_istart_64.png
   :width: 600 px
   :figwidth: 600 px
   :align: center

A subset of galaxies with *50th perentile* in brightness:

.. figure:: _static/fake_real_mosaic_istart_112.png
   :width: 600 px
   :figwidth: 600 px
   :align: center

A subset of galaxies with *25th perentile* in brightness:

.. figure:: _static/fake_real_mosaic_istart_208.png
   :width: 600 px
   :figwidth: 600 px
   :align: center

A subset of the *faintest* galaxies (1st percentile in brightness):

.. figure:: _static/fake_real_mosaic_istart_254.png
   :width: 600 px
   :figwidth: 600 px
   :align: center
