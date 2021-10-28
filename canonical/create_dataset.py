###############################################################################
'''This parser takes as input the .h5 file produced by create_datafile.py and
outputs a .h5 file with datapoints of the form (X, Y), which can be understood
by Keras models.'''
###############################################################################

import h5py
import numpy as np
import sys
import time
from utils import *
from constants import *

start_time = time.time()

assert sys.argv[1] in ['train', 'test', 'all']
assert sys.argv[2] in ['0', '1', 'all']

h5f = h5py.File(data_dir + 'datafile'
                + '_' + sys.argv[1] + '_' + sys.argv[2]
                + '.h5', 'r')

'''
each of these variable has 1.652 and 13.384 values
for test set and train set respectively
'''
SEQ = h5f['SEQ'][:]
STRAND = h5f['STRAND'][:]
TX_START = h5f['TX_START'][:]
TX_END = h5f['TX_END'][:]
JN_START = h5f['JN_START'][:]
JN_END = h5f['JN_END'][:]
h5f.close()

h5f2 = h5py.File(data_dir + 'dataset'
                + '_' + sys.argv[1] + '_' + sys.argv[2]
                + '.h5', 'w')

CHUNK_SIZE = 100

'''
This for loop set the label for each base in the gene sequence by create_datapoints function.
Dataset is created by using h5py
'''

for i in range(SEQ.shape[0]//CHUNK_SIZE):
    # Each dataset has CHUNK_SIZE genes
    
    if (i+1) == SEQ.shape[0]//CHUNK_SIZE:
        NEW_CHUNK_SIZE = CHUNK_SIZE + SEQ.shape[0]%CHUNK_SIZE
    else:
        NEW_CHUNK_SIZE = CHUNK_SIZE

    X_batch = []
    Y_batch = [[]]

    #this loop is executed number of genes times
    for j in range(NEW_CHUNK_SIZE):
        
        #idx from 0 to SEQ.shape[0]
        idx = i*CHUNK_SIZE + j

        # create_datapoints from utils.py
        X, Y = create_datapoints(SEQ[idx], STRAND[idx],
                                 TX_START[idx], TX_END[idx],
                                 JN_START[idx], JN_END[idx])

        X_batch.extend(X)
        Y_batch[0].extend(Y[0])

    X_batch = np.asarray(X_batch).astype('int8')
    Y_batch[0] = np.asarray(Y_batch[0]).astype('int8')

    #h5f2 has (SEQ.shape[0]//CHUNK_SIZE)*2 keys
    h5f2.create_dataset('X' + str(i), data=X_batch)
    h5f2.create_dataset('Y' + str(i), data=Y_batch)

h5f2.close()

print ("--- %s seconds ---" % (time.time() - start_time))

###############################################################################         
