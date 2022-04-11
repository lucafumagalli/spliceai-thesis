import numpy as np
import re
import sys
import time
import h5py
from constants import *

start_time = time.time()

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
        if (table_list[0] == "CFTR"):

            sequence_list = re.split('\n|\t|:|-', sequence_line)[:-1]
        
            chrom_table = table_list[2]
            chrom_sequence = sequence_list[0]
            transcr_start = table_list[4]
            transcr_end = table_list[5]
            paralog = table_list[1]

    
            assert chrom_table == chrom_sequence
            assert int(transcr_start) == int(sequence_list[1])+CL_max//2+1
            assert int(transcr_end) == int(sequence_list[2])-CL_max//2

                    
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
                + '_CFTR.h5', 'w')

print("NAME: ", NAME[0])
print("PARALOG: ", PARALOG)
print("CHROM: ", CHROM[0])
print("STRAND: ", STRAND[0])
print("TX_START: ", TX_START[0])
print("TX_END: ", TX_END[0])
print("JN_START: ", JN_START[0] )
print("JN_END: ", JN_END[0])


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

print(h5f['SEQ'])

h5f.close()

print ("--- %s seconds ---" % (time.time() - start_time))

###############################################################################


