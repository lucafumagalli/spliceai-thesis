###############################################################################
# This file contains code to test the SpliceAI model.
###############################################################################

import numpy as np
import sys
import time
import h5py
from keras.models import load_model
from utils import *
from constants import *
from tqdm import tqdm

assert int(sys.argv[1]) in [80, 400, 2000, 10000]
CL = int(sys.argv[1])

###############################################################################
# Load model and test data
###############################################################################

BATCH_SIZE = 6
version = [1, 2, 3, 4, 5]

model = [[] for v in range(len(version))]
'''
Load the 5 models and load them into model[]
'''
for v in range(len(version)):
    model[v] = load_model('Models/SpliceAI' + str(CL)
                          + '_c' + str(version[v]) + '.h5', compile=False)

h5f = h5py.File(data_dir + 'dataset_CFTR.h5', 'r')

'''
Get the number of keys  of the test dataset
'''
num_idx = len(h5f.keys())//2
print("h5f keys:", h5f.keys())

###############################################################################
'''
Model testing
'''
###############################################################################

start_time = time.time()

'''
The three neurons per output correspond to no splicing, 
splice acceptor (AG) and splice donor (GT) respectively.
'''
output_class_labels = ['Null', 'Acceptor', 'Donor']


for output_class in [1, 2]:

    Y_true = [[] for t in range(1)]
    Y_pred = [[] for t in range(1)]

    '''
    Loop over each subset of the test set
    '''
    for idx in tqdm(range(num_idx)):
        print("prova")

        '''
        Take data and their labels
        '''
        X = h5f['X' + str(idx)][:]
        Y = h5f['Y' + str(idx)][:]

        Xc, Yc = clip_datapoints(X, Y, CL, 1)

        Yps = [np.zeros(Yc[0].shape)]

        '''
        Prediction for each of the 5 model
        '''
        for v in range(len(version)):

            Yp = model[v].predict(Xc, batch_size=BATCH_SIZE)

            if not isinstance(Yp, list):
                Yp = [Yp]

            '''
                The final prediction is the average
                of the five predictions of the five
                different models
            '''
            Yps[0] += Yp[0]/len(version)
            
        is_expr = (Yc[0].sum(axis=(1,2)) >= 1)

        Y_true[0].extend(Yc[0][is_expr, :, output_class].flatten())
        Y_pred[0].extend(Yps[0][is_expr, :, output_class].flatten())

    print ("\n\033[1m%s:\033[0m" % (output_class_labels[output_class]))

    Y_true[0] = np.asarray(Y_true[0])
    Y_pred[0] = np.asarray(Y_pred[0])

    print_topl_statistics(Y_true[0], Y_pred[0])


h5f.close()

print ("--- %s seconds ---" % (time.time() - start_time))
print ("--------------------------------------------------------------")

###############################################################################


