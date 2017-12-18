"""
Credit: https://github.com/ageron/handson-ml
Adapted from their CNN ipynb in chapter 13
"""

import numpy as np
import os
from datetime import datetime
import tensorflow as tf

# Tensorboard
now = datetime.utcnow().strftime("%Y%m%d%H%M%S")
root_logdir = "tf_logs"
logdir = "{}/run-{}/".format(root_logdir, now)


height,width,channels = (64,64,6) 

conv_kwargs= dict(strides=1,
                  padding='SAME',
                  activation=tf.nn.relu)
pool_kwargs= dict(ksize= [1,2,2,1],
                  strides=[1,2,2,1],
                  padding='VALID')

with tf.name_scope("inputs"):
    X = tf.placeholder(tf.float32, shape=[None,height,width,channels], name="X") #training data shape
    #X_reshaped = tf.reshape(X, shape=[-1, height, width, channels])
    y = tf.placeholder(tf.int32, shape=[None], name="y") 



# 64x64
with tf.name_scope("layer1"):
    conv1 = tf.layers.conv2d(X, filters=3*channels, kernel_size=7,
                             **conv_kwargs)
    pool1 = tf.nn.avg_pool(conv1, **pool_kwargs)

# 32x32
with tf.name_scope("layer2"):
    conv2 = tf.layers.conv2d(pool1, filters=6*channels, kernel_size=7,
                             **conv_kwargs)
    pool2 = tf.nn.avg_pool(conv2, **pool_kwargs)

# 16x16
with tf.name_scope("layer3"):
    conv3 = tf.layers.conv2d(pool2, filters=9*channels, kernel_size=7,
                             **conv_kwargs)
    pool3 = tf.nn.avg_pool(conv3, **pool_kwargs)
    # next is fc
    pool3_flat = tf.reshape(pool3, 
                            shape=[-1, pool3.shape[1] * pool3.shape[2] * pool3.shape[3]])


with tf.name_scope("fc"):
    fc = tf.layers.dense(pool3_flat, 64, activation=tf.nn.relu, name="fc")

with tf.name_scope("output"):
    logits = tf.layers.dense(fc, 2, name="output") # 2 classes
    Y_proba = tf.nn.softmax(logits, name="Y_proba")

with tf.name_scope("train"):
    xentropy = tf.nn.sparse_softmax_cross_entropy_with_logits(logits=logits, labels=y)
    loss = tf.reduce_mean(xentropy)
    optimizer = tf.train.AdamOptimizer()
    training_op = optimizer.minimize(loss)

# Tensorboard
loss_summary = tf.summary.scalar('LOSS', loss)
file_writer = tf.summary.FileWriter(logdir, tf.get_default_graph())


with tf.name_scope("eval"):
    correct = tf.nn.in_top_k(logits, y, 1)
    accuracy = tf.reduce_mean(tf.cast(correct, tf.float32))

with tf.name_scope("init_and_save"):
    init = tf.global_variables_initializer()
    saver = tf.train.Saver()


# Sanity check
def get_data_dir(brick):
    if os.environ['HOME'] == '/Users/kaylan1':
        dr= os.path.join(os.environ['HOME'],'Downloads')
    else: 
        dr= os.path.join(os.environ['CSCRATCH'],'obiwan_out')
    return os.path.join(dr,'dr5_testtrain','testtrain',brick[:3],brick) 
    
dr= get_data_dir('1211p060')
Xtrain= np.load(os.path.join(dr,'xtrain_1.npy'))
Ytrain= np.load(os.path.join(dr,'ytrain_1.npy'))

batch_size=32
X_batch,y_batch= Xtrain[:batch_size,...],Ytrain[:batch_size]

with tf.Session() as sess:
    init.run()
    print(sess.run(X, feed_dict={X: X_batch}).shape)
    print(sess.run(conv1, feed_dict={X: X_batch}).shape)
    print(sess.run(pool1, feed_dict={X: X_batch}).shape)
    print(sess.run(conv2, feed_dict={X: X_batch}).shape)
    print(sess.run(pool2, feed_dict={X: X_batch}).shape)
    print(sess.run(conv3, feed_dict={X: X_batch}).shape)
    print(sess.run(pool3, feed_dict={X: X_batch}).shape)
    print(sess.run(pool3_flat, feed_dict={X: X_batch}).shape)
    print(sess.run(fc, feed_dict={X: X_batch}).shape)
    print(sess.run(logits, feed_dict={X: X_batch}).shape) 

# Train
n_epochs = 8
batch_size = 32

with tf.Session() as sess:
    init.run()
    for epoch in range(n_epochs):
        for cnt in range(Xtrain.shape[0] // batch_size):
            slc= slice(cnt*batch_size,(cnt+1)*batch_size)
            X_batch, y_batch= Xtrain[slc,...],Ytrain[slc]
            sess.run(training_op, feed_dict={X: X_batch, y: y_batch})
        acc_train = accuracy.eval(feed_dict={X: X_batch, y: y_batch})
        print(epoch, "Train accuracy:", acc_train)

        save_path = saver.save(sess, "./cnn_model")
