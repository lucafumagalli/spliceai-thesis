###############################################################################
'''
This file contains the code to train the SpliceAI model.
'''
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
import matplotlib.pyplot as plt
from tqdm import tqdm
import pandas as pd


assert int(sys.argv[1]) in [80, 400, 2000, 10000]

###############################################################################
'''
Creation of the model with the hyper-parameters depending of the argv[1], the number of nucleotides.
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

accuracy_donor_train = []
accuracy_donor_val = []
accuracy_acceptor_train = []
accuracy_acceptor_val = []
df_accuracy_donor_train = pd.DataFrame(columns=['k0.5','k1','k2','k4','AUPRC','t1','t2','t3','t4','#idx_true'])
df_accuracy_donor_val = df_accuracy_donor_train.copy()
df_accuracy_acceptor_train = df_accuracy_donor_train.copy()
df_accuracy_acceptor_val = df_accuracy_donor_train.copy()



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

'''
Take indeces for train and validation from the random permutation of dataset indeces
'''
idx_train = idx_all[:int(0.9*num_idx)]
idx_valid = idx_all[int(0.9*num_idx):]

'''
The model is trained for 10 epochs
'''
epochs = 10
EPOCH_NUM = epochs*len(idx_train)

start_time = time.time()

'''
In each epoch get a index from idx_train and use it to train the model on  h5f[idx] values.
'''
for epoch_num in tqdm(range(EPOCH_NUM)):

    idx = np.random.choice(idx_train)

    X = h5f['X' + str(idx)][:]
    Y = h5f['Y' + str(idx)][:]

    Xc, Yc = clip_datapoints(X, Y, CL, N_GPUS) 
    model_m.fit(Xc, Yc, batch_size=BATCH_SIZE, verbose=0)

    '''
    When the length of idx_train is a multiple of the number of the successive epoch,
    print metrics(see utils.py for details) for validation and training set
    '''
    if (epoch_num+1) % len(idx_train) == 0:
        
        print ("--------------------------------------------------------------")
        print ("\n\033[1mValidation set metrics:\033[0m")

        Y_true_1 = [[]]
        Y_true_2 = [[]]
        Y_pred_1 = [[]]
        Y_pred_2 = [[]]

        '''
        Predict the validation set
        '''
        for idx in idx_valid:

            X = h5f['X' + str(idx)][:]
            Y = h5f['Y' + str(idx)][:]

            Xc, Yc = clip_datapoints(X, Y, CL, N_GPUS)
            '''
            Prediction of a subset of the validation set
            '''
            Yp = model_m.predict(Xc, batch_size=BATCH_SIZE)

            if not isinstance(Yp, list):
                Yp = [Yp]


            is_expr = (Yc[0].sum(axis=(1,2)) >= 1)

            Y_true_1[0].extend(Yc[0][is_expr, :, 1].flatten())
            Y_true_2[0].extend(Yc[0][is_expr, :, 2].flatten())
            Y_pred_1[0].extend(Yp[0][is_expr, :, 1].flatten())
            Y_pred_2[0].extend(Yp[0][is_expr, :, 2].flatten())

        
        '''
        Print the statistics for the validation set
        '''
        print ("\n\033[1mAcceptor:\033[0m")
        accuracy_list, threshold_list, idx_true, auprc = print_topl_statistics(np.asarray(Y_true_1[0]),np.asarray(Y_pred_1[0]))
        accuracy_acceptor_val.append(accuracy_list[3])
        accuracy_list.append(auprc)
        accuracy_list.extend(threshold_list)
        accuracy_list.append(idx_true)
        df_accuracy_acceptor_val.loc[len(df_accuracy_acceptor_val)] = accuracy_list

        print ("\n\033[1mDonor:\033[0m")
        accuracy_list, threshold_list, idx_true, auprc = print_topl_statistics(np.asarray(Y_true_2[0]),np.asarray(Y_pred_2[0]))
        accuracy_donor_val.append(accuracy_list[3])
        accuracy_list.append(auprc)
        accuracy_list.extend(threshold_list)
        accuracy_list.append(idx_true)
        df_accuracy_donor_val.loc[len(df_accuracy_donor_val)] = accuracy_list


        print ("\n\033[1mTraining set metrics:\033[0m")

        Y_true_1 = [[]]
        Y_true_2 = [[]]
        Y_pred_1 = [[]]
        Y_pred_2 = [[]]

        '''
        Predict the first len(idx_valid) elements of the training set
        '''
        for idx in idx_train[:len(idx_valid)]:

            X = h5f['X' + str(idx)][:]
            Y = h5f['Y' + str(idx)][:]

            Xc, Yc = clip_datapoints(X, Y, CL, N_GPUS)
            Yp = model_m.predict(Xc, batch_size=BATCH_SIZE)

            if not isinstance(Yp, list):
                Yp = [Yp]

            '''
            Create a list of boolean values where is True 
            '''
            is_expr = (Yc[0].sum(axis=(1,2)) >= 1)

            Y_true_1[0].extend(Yc[0][is_expr, :, 1].flatten())
            Y_true_2[0].extend(Yc[0][is_expr, :, 2].flatten())
            Y_pred_1[0].extend(Yp[0][is_expr, :, 1].flatten())
            Y_pred_2[0].extend(Yp[0][is_expr, :, 2].flatten())

        '''
        Print the statistics for the training set
        '''
        print ("\n\033[1mAcceptor:\033[0m")
        accuracy_list, threshold_list, idx_true, auprc = print_topl_statistics(np.asarray(Y_true_1[0]),np.asarray(Y_pred_1[0]))
        accuracy_acceptor_train.append(accuracy_list[3])
        accuracy_list.append(auprc)
        accuracy_list.extend(threshold_list)
        accuracy_list.append(idx_true)
        df_accuracy_acceptor_train.loc[len(df_accuracy_acceptor_train)] = accuracy_list

        print ("\n\033[1mDonor:\033[0m")
        accuracy_list, threshold_list, idx_true, auprc = print_topl_statistics(np.asarray(Y_true_2[0]),np.asarray(Y_pred_2[0]))
        accuracy_donor_train.append(accuracy_list[3])
        accuracy_list.append(auprc)
        accuracy_list.extend(threshold_list)
        accuracy_list.append(idx_true)
        df_accuracy_donor_train.loc[len(df_accuracy_donor_train)] = accuracy_list

        print( "Learning rate: %.5f" % (kb.get_value(model_m.optimizer.lr)))
        print ("--- %s seconds ---" % (time.time() - start_time))
        start_time = time.time()

        print ("--------------------------------------------------------------")

        model.save('./Models/SpliceAI' + sys.argv[1]
                   + '_c' + sys.argv[2] + '.h5')

        '''
        The learning rate of the optimizer was set to 0.001(default) for the first 6 epochs,
        and then reduced by a factor of 2 in every subsequent epoch
        '''
        if (epoch_num+1) >= 6*len(idx_train):
            kb.set_value(model_m.optimizer.lr,
                         0.5*kb.get_value(model_m.optimizer.lr))

#saving dataframe results
df_accuracy_donor_train.to_pickle("./donor_train" +str(sys.argv[1]) + "_" + str(sys.argv[2]) + ".pkl")
df_accuracy_donor_val.to_pickle("./donor_val" + str(sys.argv[1]) + "_" +str(sys.argv[2]) + ".pkl")
df_accuracy_acceptor_train.to_pickle("./acceptor_train " + str(sys.argv[1])+"_" +str(sys.argv[2]) + ".pkl")
df_accuracy_acceptor_val.to_pickle("./acceptor_val"+ str(sys.argv[1]) + "_"+ str(sys.argv[2]) + ".pkl")

#plot for training accuracy
fig = plt.figure()
plt.plot([ep for ep in range(epochs)],[acc for acc in accuracy_acceptor_train], label='acceptor')
plt.plot([ep for ep in range(epochs)],[don for don in accuracy_donor_train], label='donor')
plt.xlabel('epoch')
plt.ylabel('topk-accuracy')
plt.legend()
plt.xticks([x for x in range(epochs)])
fig.savefig('accuracy_train' + str(sys.argv[1]) + '_' + str(sys.argv[2]) + '.png')

#plot for validation accuracy
fig1 = plt.figure()
plt.plot([ep for ep in range(epochs)],[acc for acc in accuracy_acceptor_val], label='acceptor')
plt.plot([ep for ep in range(epochs)],[don for don in accuracy_donor_val], label='donor')
plt.xlabel('epoch')
plt.ylabel('topk-accuracy')
plt.legend()
plt.xticks([x for x in range(epochs)])
fig1.savefig('accuracy_val' + str(sys.argv[1])+ '_' + str(sys.argv[2]) + '.png')


#plot for train/val acceptor
fig2 = plt.figure()
plt.plot([ep for ep in range(epochs)],[acc for acc in accuracy_acceptor_val], label='validation')
plt.plot([ep for ep in range(epochs)],[don for don in accuracy_acceptor_train], label='training')
plt.xlabel('epoch')
plt.ylabel('topk-accuracy')
plt.legend()
plt.xticks([x for x in range(epochs)])
fig2.savefig('acceptor_train_val' + str(sys.argv[1])+ '_'+ str(sys.argv[2]) + '.png')

#plot for train/val donor
fig3 = plt.figure()
plt.plot([ep for ep in range(epochs)],[acc for acc in accuracy_donor_val], label='validation')
plt.plot([ep for ep in range(epochs)],[don for don in accuracy_donor_train], label='training')
plt.xlabel('epoch')
plt.ylabel('topk-accuracy')
plt.legend()
plt.xticks([x for x in range(epochs)])
fig3.savefig('donor_train_val' + str(sys.argv[1]) + '_ ' + str(sys.argv[2]) + '.png')

h5f.close()
        
###############################################################################
