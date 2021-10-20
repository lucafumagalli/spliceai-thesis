# SpliceAI

Download ```hg19.fa``` file (hg19/GRCh37) from [GRCh37/hg19](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz) and save it in ```canonical``` folder.

## Datasets creation
Use a python2 env and install requirements:   
```
pip install -r requirements/python2.txt
```
From ```canonical``` folder run the following command to create the sequence file:

```
./grabsequence.sh
```
Then, to create datasets run:

```
./create_datafile.sh  
./create_dataset.sh
``` 

