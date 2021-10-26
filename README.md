# SpliceAI

Download ```hg19.fa``` file (hg19/GRCh37) from [GRCh37/hg19](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz) and save it in ```canonical``` folder.

## Requirements
Use a python2 env and install requirements:
```
pip install -r requirements.txt
```

## Datasets creation
From ```canonical``` folder run the following command to create the sequence file:

```
./grabsequence.sh
```
Then, to create datasets run:

```
./create_datafile.sh  
./create_dataset.sh
``` 
## Train models
In ```canonical``` folder create a ```Outputs``` folder for saving outputs of different training, create also a ```Models``` folder for saving the models.

For training the 5 models with 10.000 nucleotides run the following commands from the ```canonical``` folder:
```
./script_train.sh 10000 1
./script_train.sh 10000 2
./script_train.sh 10000 3
./script_train.sh 10000 4
./script_train.sh 10000 5
```
## Test models
For training the model with 10.000 nucleotides, from the ```canonical``` folder run:
```
./script_test.sh 10000
```