"""
Trains a CNN on fake and real galaxy images using TensorFlow.

Adapted from https://github.com/ageron/handson-ml
"""

import numpy as np
import os
from datetime import datetime
from glob import glob

import tensorflow as tf

def get_indir(nersc=False):
    '''Returns path to dr5_testtrain directory'''
    if nersc:
        return os.path.join('/global/cscratch1/sd/kaylanb','obiwan_out')
    else:
        return os.path.join(os.environ['HOME'],'Downloads')

def get_outdir(nersc=False,knl=False):
    '''Where to write ckpt and log files'''
    if (nersc) & (knl):
        return os.path.join('/global/cscratch1/sd/kaylanb','obiwan_out','cnn_knl')
    elif (nersc) & (not knl):
        return os.path.join('/global/cscratch1/sd/kaylanb','obiwan_out','cnn')
    else:
        return os.path.join(os.environ['HOME'],'Downloads','cnn')

def get_xtrain_fns(brick,indir):
    search= os.path.join(indir,'dr5_testtrain','testtrain',
                         brick[:3],brick,'xtrain_*.npy')
    xtrain_fns= glob(search)
    if len(xtrain_fns) == 0:
        raise IOError('No training data found matching: %s' % search)
    return xtrain_fns

#def BatchGen(X,y,brick,batch_size=32):
def BatchGen(brick,indir,batch_size=32):
    fns= get_xtrain_fns(brick,indir)
    for fn in fns:
        print('Loading %s' % fn)
        X= np.load(fns[0])
        y= np.load(fns[0].replace('xtrain_','ytrain_'))
        N= X.shape[0]
        ind= np.array_split(np.arange(N),N // batch_size + 1)
        for i in ind:
            yield X[i,...],y[i].astype(np.int32) #.reshape(-1,1).astype(np.int32)

def get_bricks(fn='cnn_bricks.txt'):
    fn= os.path.join(os.path.dirname(__file__),
                     '../../../etc',fn)
    if not os.path.exists(fn):
	    raise IOError('Need to create brick list: %s' % fn)
    bricks= np.loadtxt(fn,dtype=str)
    if len(bricks.shape) == 0:
        # single line
        with open(fn,'r') as f:
            bricks= np.array([f.read().strip()])
    return bricks


def get_checkpoint(epoch,brick,outdir):
    return os.path.join(outdir,'ckpts',
                        'epoch_%s_brick_%s.ckpt' % (epoch,brick))

def bookmark_fn(outdir):
	"""Single line text file storing the epoch,brick,batch number of last ckpt"""
	return os.path.join(outdir,'ckpts',
                        'last_epoch_brick_batch.txt')

def get_bookmark(outdir):
    with open(bookmark_fn(outdir),'r') as f:
        epoch,brick,ith_batch= f.read().strip().split(' ')
    return epoch,brick,ith_batch

def get_logdir(outdir):
    now = datetime.utcnow().strftime("%Y%m%d%H%M%S")
    logdir= os.path.join(outdir,'logs')
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


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser= ArgumentParser()
    parser.add_argument('--outdir', type=str, default=None, help='optional output directory for the checkpoint and log files')
    args= parser.parse_args()

    knl=False
    config=None
    if 'isKNL' in os.environ.keys():
        # Set in slurm_job_knl.sh
        knl=True
        config= tf.ConfigProto()
        config.intra_op_parallelism_threads=os.environ['OMP_NUM_THREADS']
        assert(os.environ['OMP_NUM_THREADS'] == 68)
        config.inter_op_parallelism_threads=1

    nersc=False
    if 'CSCRATCH' in os.environ.keys():
        nersc=True
    indir= get_indir(nersc=nersc)
    outdir= get_outdir(nersc=nersc,knl=knl)
    if not args.outdir is None:
        outdir= args.outdir

    # Train
    n_epochs = 4
    batch_size = 16
    bricks= get_bricks()
    file_writer = tf.summary.FileWriter(get_logdir(outdir),
                                        tf.get_default_graph())

    first_epoch,first_brick,first_batch= '0',bricks[0],'0'
    fn= get_checkpoint(first_epoch,first_brick,outdir)+'.meta'
    if os.path.exists(fn):
        last_epoch,last_brick,last_batch= get_bookmark(outdir)
        ckpt_fn= get_checkpoint(last_epoch,last_brick,outdir)
    else:
        last_epoch,last_brick,last_batch= first_epoch,first_brick,first_batch
        ckpt_fn= None
    last_ibrick= np.where(bricks == last_brick)[0][0] #+ 1 creates bug where last break skips all epochs
    #bricks= ['1211p060']

    with tf.Session(config=config) as sess:
        if ckpt_fn is None:
            sess.run(init)
            print('Starting from scratch')
        else:
            saver.restore(sess, ckpt_fn)
            print('Restored ckpt %s' % ckpt_fn)

        batch_index= int(last_batch)
        for epoch in range(int(last_epoch),n_epochs+1):
            for ibrick,brick in enumerate(bricks):
                # Don't repeat bricks when restart from ckpt
                if ibrick < last_ibrick:
                    print('skipping: epoch,ibrick,last_ibrick', epoch,ibrick,last_ibrick)
                    continue
                data_gen= BatchGen(brick,indir,batch_size)
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
                # Save progress
                fn= get_checkpoint(epoch,brick,outdir)
                save_path = saver.save(sess, fn)
                print('Wrote ckpt %s' % fn)
                with open(bookmark_fn(outdir),'w') as f:
                    f.write('%d %s %d' % (epoch,brick,batch_index))
                print('Updated %s' % bookmark_fn(outdir))
                # Reset last_ibrick so use all bricks in next epoch
                last_ibrick= 0
