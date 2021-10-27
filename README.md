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

## Files description 
### ```constants.py```
This file sets the maximum nucleotide context length and the sequence length of the SpliceAI models.  
It sets also the path of the canonical dataset and the path of the FASTA file of the genome.
### ```utils.py```
In this file are present the functions for creating the datasets: creating datapoints, reformatting data and one hot encoding sequence.   
There is also a function for printing the statistics after the model testing.
### ```create_datafile.py```
This parser takes as input the text files canonical_dataset.txt and canonical_sequence.txt, and produces a .h5 file datafile.
### ```create_dataset.py```
This parser takes as input the .h5 file produced by create_datafile.py and outputs a .h5 file with datapoints of the form (X, Y), which can be understood by Keras models.
### ```spliceai.py```
This file has the functions to create the spliceAI model, that is represent in the figure below, depending on the number of nucleotides used a different model is created.
![alt text](https://raw.githubusercontent.com/lucafumagalli/spliceai-thesis/main/images/architectures.png?token=ACTSO2SZZRLJ7X57S7FGR3LBPF3TK)
### ```train_model.py```
This file contains the code to train the SpliceAI model.
### ```test_model.py```
Contains code to test the SpliceAI model.
### ```multi_gpu.py```


