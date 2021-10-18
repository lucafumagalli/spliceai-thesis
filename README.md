# SpliceAI

## Datasets creation
Use a python2 env and install requirements:   
```
pip install -r requirements/python2.txt
```
Use the following commands from ```canonical``` folder.   
To download genome file(hg19/GRCh37) from [UCSC](https://hgdownload.soe.ucsc.edu/downloads.html) and create sequence file run:

```
./grabsequence.sh
```
Then, to create datasets run:

```
./create_datafile.sh  
./create_dataset.sh
``` 
```utils.py``` has functions which process the information in the .h5 files and convert them into a format usable by Keras.   

```spliceai.py``` has the functions necessary to create the SpliceAI model.   

```test_model.py``` contains code to test the SpliceAI model.  

The code was tested using ```keras==2.0.5``` and ```tensorflow==1.11.0```
