#!/bin/bash

source constants.py

CLr=$((CL_max/2))
CLl=$(($CLr+1))
# First nucleotide not included by BEDtools

#from canonical_dataset.txt get the interval(first and last position) for each gene
#substracting and summing the context length and save it in temp.bed
cat $splice_table | awk -v CLl=$CLl -v CLr=$CLr '{print $3"\t"($5-CLl)"\t"($6+CLr)}' > temp.bed

#extracts sequences from genome.fa file for each of the intervals defined in temp.bed and save them in canonical_sequence.txt
bedtools getfasta -bed temp.bed -fi $ref_genome -fo $sequence -tab

rm temp.bed
