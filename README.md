# SpliceAI

## Datasets creation
Use a python2 env and install requirements:   
```
pip install -r requirements/python2.txt
```
To download genome file(hg19/GRCh37) from [UCSC](https://hgdownload.soe.ucsc.edu/downloads.html) and create sequence file run:

```
./canonical/grabsequence.sh
```
Then, use the following commands to create datasets:

```
./canonical/create_datafile.sh  
./canonical/create_dataset.sh

``` 
```utils.py``` has functions which process the information in the .h5 files and convert them into a format usable by Keras.   

```spliceai.py``` has the functions necessary to create the SpliceAI model.   

```test_model.py``` contains code to test the SpliceAI model.  

The code was tested using ```keras==2.0.5``` and ```tensorflow==1.11.0```
