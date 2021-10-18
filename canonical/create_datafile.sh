#!/bin/bash

#create datafile for training and testing

python2.7 create_datafile.py train all
python2.7 create_datafile.py test 0