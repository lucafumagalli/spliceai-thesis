###############################################################################
# This file contains the code to train the SpliceAI model.
###############################################################################

import numpy as np
import sys
import time
import h5py
import keras.backend as kb
import tensorflow as tf
from spliceai import *
from utils import *
from multi_gpu import *
from constants import * 

assert int(sys.argv[1]) in [80, 400, 2000, 10000]

###############################################################################
'''
Creation of the model with the right hyper-parameters depending of the argv[1], the number of nucleotides.
'''
###############################################################################

L = 32
N_GPUS = 2

# Hyper-parameters:
# L: Number of convolution kernels
# W: Convolution window size in each residual unit
# AR: Atrous rate in each residual unit

#case with 80 nt
if int(sys.argv[1]) == 80:
    W = np.asarray([11, 11, 11, 11])
    AR = np.asarray([1, 1, 1, 1])
    BATCH_SIZE = 18*N_GPUS
#case with 400 nt
elif int(sys.argv[1]) == 400:
    W = np.asarray([11, 11, 11, 11, 11, 11, 11, 11])
    AR = np.asarray([1, 1, 1, 1, 4, 4, 4, 4])
    BATCH_SIZE = 18*N_GPUS
#cae with 2.000 nt
elif int(sys.argv[1]) == 2000:
    W = np.asarray([11, 11, 11, 11, 11, 11, 11, 11,
                    21, 21, 21, 21])
    AR = np.asarray([1, 1, 1, 1, 4, 4, 4, 4,
                     10, 10, 10, 10])
    BATCH_SIZE = 12*N_GPUS
#case with 10.000 nt
elif int(sys.argv[1]) == 10000:
    W = np.asarray([11, 11, 11, 11, 11, 11, 11, 11,
                    21, 21, 21, 21, 41, 41, 41, 41])
    AR = np.asarray([1, 1, 1, 1, 4, 4, 4, 4,
                     10, 10, 10, 10, 25, 25, 25, 25])
    BATCH_SIZE = 6*N_GPUS

#the context nucleotides
CL = 2 * np.sum(AR*(W-1))
assert CL <= CL_max and CL == int(sys.argv[1])
print ("\033[1mContext nucleotides: %d\033[0m" % (CL))
print ("\033[1mSequence length (output): %d\033[0m" % (SL))

#model creation
model = SpliceAI(L, W, AR)
model_m = make_parallel(model, N_GPUS)
#loss function defined in spliceai.py
model_m.compile(loss=categorical_crossentropy_2d, optimizer='adam')

###############################################################################
# Training and validation
###############################################################################

h5f = h5py.File(data_dir + 'dataset' + '_' + 'train'
                + '_' + 'all' + '.h5', 'r')

num_idx = len(h5f.keys())//2
idx_all = np.random.permutation(num_idx)
#take indeces for train and test from the random permutation of dataset indeces
idx_train = idx_all[:int(0.9*num_idx)]
idx_valid = idx_all[int(0.9*num_idx):]

EPOCH_NUM = 10*len(idx_train)

start_time = time.time()


for epoch_num in range(EPOCH_NUM):

    idx = np.random.choice(idx_train)

    X = h5f['X' + str(idx)][:]
    Y = h5f['Y' + str(idx)][:]

    Xc, Yc = clip_datapoints(X, Y, CL, N_GPUS) 
    model_m.fit(Xc, Yc, batch_size=BATCH_SIZE, verbose=0)

    #print metrics of validation for check performances
    if (epoch_num+1) % len(idx_train) == 0:
        # Printing metrics (see utils.py for details)

        print ("--------------------------------------------------------------")
        print ("\n\033[1mValidation set metrics:\033[0m")

        Y_true_1 = [[]]
        Y_true_2 = [[]]
        Y_pred_1 = [[]]
        Y_pred_2 = [[]]

        #predict value for validation indices
        for idx in idx_valid:

            X = h5f['X' + str(idx)][:]
            Y = h5f['Y' + str(idx)][:]

            Xc, Yc = clip_datapoints(X, Y, CL, N_GPUS)
            Yp = model_m.predict(Xc, batch_size=BATCH_SIZE)

            if not isinstance(Yp, list):
                Yp = [Yp]


            is_expr = (Yc[0].sum(axis=(1,2)) >= 1)

            Y_true_1[0].extend(Yc[0][is_expr, :, 1].flatten())
            Y_true_2[0].extend(Yc[0][is_expr, :, 2].flatten())
            Y_pred_1[0].extend(Yp[0][is_expr, :, 1].flatten())
            Y_pred_2[0].extend(Yp[0][is_expr, :, 2].flatten())

        print ("\n\033[1mAcceptor:\033[0m")
        print_topl_statistics(np.asarray(Y_true_1[0]),
                                np.asarray(Y_pred_1[0]))

        print ("\n\033[1mDonor:\033[0m")
        print_topl_statistics(np.asarray(Y_true_2[0]),
                                  np.asarray(Y_pred_2[0]))

        print ("\n\033[1mTraining set metrics:\033[0m")

        Y_true_1 = [[]]
        Y_true_2 = [[]]
        Y_pred_1 = [[]]
        Y_pred_2 = [[]]

        for idx in idx_train[:len(idx_valid)]:

            X = h5f['X' + str(idx)][:]
            Y = h5f['Y' + str(idx)][:]

            Xc, Yc = clip_datapoints(X, Y, CL, N_GPUS)
            Yp = model_m.predict(Xc, batch_size=BATCH_SIZE)

            if not isinstance(Yp, list):
                Yp = [Yp]

            is_expr = (Yc[0].sum(axis=(1,2)) >= 1)

            Y_true_1[0].extend(Yc[0][is_expr, :, 1].flatten())
            Y_true_2[0].extend(Yc[0][is_expr, :, 2].flatten())
            Y_pred_1[0].extend(Yp[0][is_expr, :, 1].flatten())
            Y_pred_2[0].extend(Yp[0][is_expr, :, 2].flatten())

        print ("\n\033[1mAcceptor:\033[0m")
        print_topl_statistics(np.asarray(Y_true_1[0]),
                                np.asarray(Y_pred_1[0]))

        print ("\n\033[1mDonor:\033[0m")
        print_topl_statistics(np.asarray(Y_true_2[0]),
                                np.asarray(Y_pred_2[0]))

        print( "Learning rate: %.5f" % (kb.get_value(model_m.optimizer.lr)))
        print ("--- %s seconds ---" % (time.time() - start_time))
        start_time = time.time()

        print ("--------------------------------------------------------------")

        model.save('./Models/SpliceAI' + sys.argv[1]
                   + '_c' + sys.argv[2] + '.h5')

        if (epoch_num+1) >= 6*len(idx_train):
            kb.set_value(model_m.optimizer.lr,
                         0.5*kb.get_value(model_m.optimizer.lr))
            # Learning rate decay

h5f.close()
        
###############################################################################

