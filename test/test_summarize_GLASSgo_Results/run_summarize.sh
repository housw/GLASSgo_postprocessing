#!/bin/bash


source activate py3k

# convert GLASSgo fasta to tsv
./GLASSgo2tsv.py 01_GLASSgo_Results/GLASSgo_output_HG001_01345__HG001_03252.fa --force

# summarize count distribution, ignore leaves that have no rank, mainly strains
#./summarize_GLASSgo_sRNA_distribution.py 01_GLASSgo_Results/ --force --rank

# summarize features
#./summarize_GLASSgo_sRNA_features.py 01_GLASSgo_Results/ --force

source deactivate
