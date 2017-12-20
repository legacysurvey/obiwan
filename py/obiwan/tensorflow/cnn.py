"""
Trains a CNN on 64x64x6 image cutouts of real and fake galaxies

Credit: Adapted from Chp 13 of https://github.com/ageron/handson-ml
"""

import numpy as np
import os
from datetime import datetime
from glob import glob

import tensorflow as tf

def get_xtrain_fns(brick):
    if os.environ['HOME'] == '/Users/kaylan1':
        dr= os.path.join(os.environ['HOME'],'Downloads')
    else: 
        dr= os.path.join(os.environ['CSCRATCH'],'obiwan_out')
    xtrain_fns= glob( os.path.join(dr,'dr5_testtrain','testtrain',
                                   brick[:3],brick,'xtrain_*.npy'))
    assert(len(xtrain_fns) > 0)
    return xtrain_fns

#def BatchGen(X,y,brick,batch_size=32):
def BatchGen(brick,batch_size=32):
    fns= get_xtrain_fns(brick)
    for fn in fns:
        print('Loading %s' % fn)
        X= np.load(fns[0])
        y= np.load(fns[0].replace('xtrain_','ytrain_'))
        N= X.shape[0]
        ind= np.array_split(np.arange(N),N // batch_size + 1)
        for i in ind:
            yield X[i,...],y[i].astype(np.int32) #.reshape(-1,1).astype(np.int32)

def get_bricks():
    bricks= np.loadtxt('cnn_bricks.txt',dtype=str)
    if len(bricks.shape) == 0:
        # single line
        with open('cnn_bricks.txt','r') as f:
            bricks= [f.read().strip()]
    return bricks


def get_checkpoint(epoch,brick):
    dr= os.path.join(os.environ['HOME'],'Downloads','cnn','ckpts')
    return os.path.join(dr,'epoch_%s_brick_%s.ckpt' % (epoch,brick))

def last_epoch_and_brick_fn():
    dr= os.path.join(os.environ['HOME'],'Downloads','cnn','ckpts')
    fn= os.path.join(dr,'last_epoch_brick.txt')
    return fn

def last_epoch_and_brick():
    dr= os.path.join(os.environ['HOME'],'Downloads','cnn','ckpts')
    with open(last_epoch_and_brick_fn(),'r') as f:
        epoch,brick= f.read().strip().split(' ')
    return epoch,brick
        
def first_epoch_and_brick():
    brick= np.loadtxt('cnn_bricks.txt',dtype=str)[0]
    return '0',brick

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






# Train
n_epochs = 4
batch_size = 16
bricks= get_bricks()
file_writer = tf.summary.FileWriter(get_logdir(), 
                                    tf.get_default_graph())

first_epoch,first_brick= '0',bricks[0]
fn= get_checkpoint(first_epoch,first_brick)+'.meta'
if os.path.exists(fn):
    last_epoch,last_brick= last_epoch_and_brick()
    ckpt_fn= get_checkpoint(last_epoch,last_brick)
else:
    last_epoch,last_brick= first_epoch,first_brick
    ckpt_fn= None

#bricks= bricks[np.where(bricks == last_brick)[0]+1:]
#raise ValueError
#bricks= ['1211p060']

with tf.Session() as sess:
    if ckpt_fn is None: 
        sess.run(init)
        print('Starting from scratch')
    else:
        saver.restore(sess, ckpt_fn)
        print('Restored ckpt %s' % ckpt_fn)
    batch_index=0
    for epoch in range(int(last_epoch),n_epochs+1):
        for brick in bricks:
            data_gen= BatchGen(brick,batch_size)
            for X_,y_ in data_gen:
                sess.run(training_op, feed_dict={X: X_, y: y_})
                batch_index+=1
                if batch_index % 2 == 0:
                    step = batch_index
                    file_writer.add_summary(loss_summary.eval(feed_dict={X: X_, y: y_}), 
                                            step)
                    file_writer.add_summary(accur_summary.eval(feed_dict={X: X_, y: y_}), 
                                            step)
                    
            acc_train = accuracy.eval(feed_dict={X: X_, y: y_})
            print(epoch, "Train accuracy:", acc_train)
            fn= get_checkpoint(epoch,brick)
            save_path = saver.save(sess, fn)
            print('Wrote ckpt %s' % fn)
            #if os.path.exists(last_epoch_and_brick_fn()):
            #    os.remove(last_epoch_and_brick_fn())
            with open(last_epoch_and_brick_fn(),'w') as f:
                f.write('%d %s' % (epoch,brick))
            print('Updated %s' % last_epoch_and_brick_fn())

