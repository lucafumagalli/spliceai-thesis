###############################################################################
'''This parser takes as input the text files canonical_dataset.txt and 
canonical_sequence.txt, and produces a .h5 file datafile_{}_{}.h5,
which will be later processed to create dataset_{}_{}.h5. The file
dataset_{}_{}.h5 will have datapoints of the form (X,Y), and can be
understood by Keras models.'''
###############################################################################

import numpy as np
import re
import sys
import time
import h5py
from constants import *

start_time = time.time()

'''
argv[1] identifies train or test set for taking the corresponding chromosomes
argv[2] identifies if paralogs need to be in the dataset
'''

assert sys.argv[1] in ['train', 'test', 'all']
assert sys.argv[2] in ['0', '1', 'all']

if sys.argv[1] == 'train':
    CHROM_GROUP = ['chr11', 'chr13', 'chr15', 'chr17', 'chr19', 'chr21',
                   'chr2', 'chr4', 'chr6', 'chr8', 'chr10', 'chr12',
                   'chr14', 'chr16', 'chr18', 'chr20', 'chr22', 'chrX', 'chrY']
elif sys.argv[1] == 'test':
    CHROM_GROUP = ['chr1', 'chr3', 'chr5', 'chr7', 'chr9']
else:
    CHROM_GROUP = ['chr1', 'chr3', 'chr5', 'chr7', 'chr9',
                   'chr11', 'chr13', 'chr15', 'chr17', 'chr19', 'chr21',
                   'chr2', 'chr4', 'chr6', 'chr8', 'chr10', 'chr12',
                   'chr14', 'chr16', 'chr18', 'chr20', 'chr22', 'chrX', 'chrY']

###############################################################################

NAME = []      # Gene symbol
PARALOG = []   # 0 if no paralogs exist, 1 otherwise
CHROM = []     # Chromosome number
STRAND = []    # Strand in which the gene lies (+ or -)
TX_START = []  # Position where transcription starts
TX_END = []    # Position where transcription ends
JN_START = []  # Positions where canonical exons end
JN_END = []    # Positions where canonical exons start
SEQ = []       # Nucleotide sequence

#open the sequence file
dna_sequence = open(sequence, 'r')

#open the annotations file
with open(splice_table, 'r') as annotations_table:
    for table_line in annotations_table:

        sequence_line = dna_sequence.readline()

        table_list = re.split('\n|\t', table_line)[:-1]
        sequence_list = re.split('\n|\t|:|-', sequence_line)[:-1]
    
        chrom_table = table_list[2]
        chrom_sequence = sequence_list[0]
        transcr_start = table_list[4]
        transcr_end = table_list[5]
        paralog = table_list[1]

 
        assert chrom_table == chrom_sequence
        assert int(transcr_start) == int(sequence_list[1])+CL_max//2+1
        assert int(transcr_end) == int(sequence_list[2])-CL_max//2

        # check if the chromosome need to be in the dataset
        if (chrom_table not in CHROM_GROUP):
            continue
        
        # check if the gene is a paralog
        if (sys.argv[2] != paralog) and (sys.argv[2] != 'all'):
            continue
        
        
        NAME.append(table_list[0])
        PARALOG.append(int(paralog))
        CHROM.append(chrom_table)
        STRAND.append(table_list[3])
        TX_START.append(transcr_start)
        TX_END.append(transcr_end)
        JN_START.append(table_list[6::2])
        JN_END.append(table_list[7::2])
        SEQ.append(sequence_list[3])

dna_sequence.close()
annotations_table.close()

###############################################################################

h5f = h5py.File(data_dir + 'datafile' 
                + '_' + sys.argv[1] + '_' + sys.argv[2]
                + '.h5', 'w')

#each of these has len: 1.652 for testing and len: 13.384 for training
h5f.create_dataset('NAME', data=np.asarray(NAME))
h5f.create_dataset('PARALOG', data=np.asarray(PARALOG))
h5f.create_dataset('CHROM', data=np.asarray(CHROM))
h5f.create_dataset('STRAND', data=np.asarray(STRAND))
h5f.create_dataset('TX_START', data=np.asarray(TX_START))
h5f.create_dataset('TX_END', data=np.asarray(TX_END))
h5f.create_dataset('JN_START', data=np.asarray(JN_START))
h5f.create_dataset('JN_END', data=np.asarray(JN_END))
h5f.create_dataset('SEQ', data=np.asarray(SEQ)) 

h5f.close()

print ("--- %s seconds ---" % (time.time() - start_time))

###############################################################################

