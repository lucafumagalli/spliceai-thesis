#!/bin/bash

#create datafile for training and testing
if [ ! -f hg19.fa ]; then
    wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
    gzip -d hg19.fa.gz
fi

python2.7 create_datafile.py train all
python2.7 create_datafile.py test 0