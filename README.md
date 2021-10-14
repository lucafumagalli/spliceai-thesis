# SpliceAI
Download the ```genome.fa``` file (hg19/GRCh37) from [GRCh37/hg19](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz) in ```/canonical``` folder

Then, use the following commands:

```
./grab_sequence.sh

python create_datafile.py train all
python create_datafile.py test 0

python create_dataset.py train all
python create_dataset.py test 0

``` 
```utils.py``` has functions which process the information in the .h5 files and convert them into a format usable by Keras.   

```spliceai.py``` has the functions necessary to create the SpliceAI model.   

```test_model.py``` contains code to test the SpliceAI model.  

The code was tested using ```keras==2.0.5``` and ```tensorflow==1.11.0```
