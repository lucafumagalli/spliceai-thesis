###############################################################################
'''
This file has the functions necessary to create the SpliceAI model.
The structure of the model depends on the number of nucleotides of input.
The basic unit of the SpliceAI architectures is a residual block,
which consists of batch-normalization layers,ReLU
and convolutional units organized in a specific manner.

The output of the models consists of three scores which sum to one,
corresponding to the probability of the position of interest being a splice acceptor, splice donor, and neither.

For creating this model is not used a sequential model, but functional API instead.
'''
###############################################################################

from keras.models import Model
from keras.layers import Input
from keras.layers.core import Activation
from keras.layers.convolutional import Conv1D, Cropping1D
from keras.layers.normalization import BatchNormalization
from keras.layers.merge import add
import keras.backend as kb
import numpy as np


def ResidualUnit(l, w, ar):
    # Residual unit proposed in "Identity mappings in Deep Residual Networks"
    # by He et al.

    #conv layers use a dilation rate to increment the receptive field
    def f(input_node):
        bn1 = BatchNormalization()(input_node)
        act1 = Activation('relu')(bn1)
        conv1 = Conv1D(l, w, dilation_rate=ar, padding='same')(act1)
        bn2 = BatchNormalization()(conv1)
        act2 = Activation('relu')(bn2)
        conv2 = Conv1D(l, w, dilation_rate=ar, padding='same')(act2)
        output_node = add([conv2, input_node])

        return output_node

    return f


def SpliceAI(L, W, AR):
    # L: Number of convolution kernels
    # W: Convolution window size in each residual unit
    # AR: Atrous rate in each residual unit

    assert len(W) == len(AR)

    CL = 2 * np.sum(AR*(W-1))

    input0 = Input(shape=(None, 4))
    conv = Conv1D(L, 1)(input0)
    skip = Conv1D(L, 1)(conv)

    #create a residual unit for each value of the hyperparameters in train_model.py
    for i in range(len(W)):
        ws = (W[i],)
        ar = (AR[i], )
        conv = ResidualUnit(L, ws, ar)(conv)
        
        #every 4 residual unit create a skip connection to the output
        if (((i+1) % 4 == 0) or ((i+1) == len(W))):
            dense = Conv1D(L, 1)(conv)
            #Layer that adds a list of inputs.
            skip = add([skip, dense])

    skip = Cropping1D(int(CL/2))(skip)

    output0 = [[]]
    output0[0] = Conv1D(3, 1, activation='softmax')(skip)
    
    model = Model(inputs=input0, outputs=output0)
    return model


def categorical_crossentropy_2d(y_true, y_pred):
    # Standard categorical cross entropy for sequence outputs

    return - kb.mean(y_true[:, :, 0]*kb.log(y_pred[:, :, 0]+1e-10)
                   + y_true[:, :, 1]*kb.log(y_pred[:, :, 1]+1e-10)
                   + y_true[:, :, 2]*kb.log(y_pred[:, :, 2]+1e-10))