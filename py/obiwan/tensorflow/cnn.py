"""
Credit: https://github.com/ageron/handson-ml
Adapted from the Chp 13 ipynb
"""

import numpy as np
import os
from datetime import datetime
import tensorflow as tf

def BatchGen(X,y,batch_size=32):
    # if not perfect divide, will drop extra training instances
    N= X.shape[0]
    ind= np.array_split(np.arange(N),N // batch_size)
    for i in ind:
        yield X[i,...],y[i].astype(np.int32) #.reshape(-1,1).astype(np.int32)

def get_xtrain_fns(brick):
    if os.environ['HOME'] == '/Users/kaylan1':
        dr= os.path.join(os.environ['HOME'],'Downloads')
    else: 
        dr= os.path.join(os.environ['CSCRATCH'],'obiwan_out')
    xtrain_fns= glob( os.path.join(dr,'dr5_testtrain','testtrain',
                                   brick[:3],brick,'xtrain_*.npy'))
    if len(xtrain_fns) > 0:
        return xtrain_fns
    else:
        return None


def get_checkpoint_fns():
    dr= os.path.join(os.environ['HOME'],'Downloads','cnn') 
    return os.path.join(dr,'check.ckpt'),
           os.path.join(dr,'final.ckpt')

def get_logdir():
    now = datetime.utcnow().strftime("%Y%m%d%H%M%S")
    logdir= os.path.join(os.environ['HOME'],'Downloads',
                         "cnn",'logs')
    return os.path.join(logdir,"{}/run-{}/".format(logdir, now))



height,width,channels = (64,64,6) 

conv_kwargs= dict(strides=1,
                  padding='SAME',
                  activation=tf.nn.relu)
pool_kwargs= dict(ksize= [1,2,2,1],
                  strides=[1,2,2,1],
                  padding='VALID')

with tf.name_scope("inputs"):
    X = tf.placeholder(tf.float32, shape=[None,height,width,channels], name="X")
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

with tf.name_scope("eval"):
    correct = tf.nn.in_top_k(logits, y, 1)
    accuracy = tf.reduce_mean(tf.cast(correct, tf.float32))

with tf.name_scope("init_and_save"):
    init = tf.global_variables_initializer()
    saver = tf.train.Saver()

init = tf.global_variables_initializer()
saver = tf.train.Saver()

loss_summary= tf.summary.scalar('loss', loss)
accur_summary = tf.summary.scalar('accuracy', accuracy)




#bricks= np.loadtxt('../../../etc/bricks_dr5_grz.txt',
#                   dtype=str)
xtrain_fns= get_xtrain_fns('1211p060')
assert(not xtrain_fns is None)
xtrain= np.load(xtrain_fns[0])
ytrain= np.load(xtrain_fns[0].replace('xtrain_','ytrain_'))

# Train
n_epochs = 4
batch_size = 32
n_batches= ytrain.shape[0]//batch_size + 1
file_writer = tf.summary.FileWriter(get_logdir(), 
                                    tf.get_default_graph())
check_fn,final_fn= get_checkpoint_fns()

with tf.Session() as sess:
    if os.path.exists(final_fn):
        saver.restore(sess, final_fn)
    else: 
        sess.run(init)
    for epoch in range(n_epochs):
        data_gen= BatchGen(xtrain,ytrain,batch_size)
        batch_index=0
        for X_,y_ in data_gen:
            sess.run(training_op, feed_dict={X: X_, y: y_})
            batch_index+=1
            if i % 1 == 0:
                step = epoch * n_batches + batch_index
                file_writer.add_summary(loss_summary.eval(feed_dict={X: X_, y: y_}), 
                                        step)
                file_writer.add_summary(accur_summary.eval(feed_dict={X: X_, y: y_}), 
                                        step)
                
        acc_train = accuracy.eval(feed_dict={X: X_, y: y_})
        print(epoch, "Train accuracy:", acc_train)
        if epoch % 2 == 0:
            save_path = saver.save(sess, check_fn)
    save_path = saver.save(sess, final_fn)

