# SpliceAI 
First, update the ```constants.py``` file:
ref_genome: path of the ```genome.fa``` file (hg19/GRCh37), can be downloaded from [GRCh37/hg19](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz)   

Then, use the following commands:

```
./grab_sequence.sh

python create_datafile.py train all
python create_datafile.py test 0

python create_dataset.py train all
python create_dataset.py test 0

```
```create_dataset.py``` and ```create_datafile.py``` need to be tested with python 2.7.
