#!/bin/bash

# TODO: something wrong with the Bacillales order, see below 
# 1385    Bacillales      order   0       0       312

# summarize count distribution, ignore leaves that have no rank, mainly strains
./summarize_GLASSgo_sRNA_distribution.py 01_GLASSgo_Results/ --force --rank

# summarize features
#./summarize_GLASSgo_sRNA_features.py 01_GLASSgo_Results/ --force
