# -*- coding: utf-8 -*-
"""
Created on Sat Aug 14 07:47:03 2021

@author: jeff
"""

import pandas as pd
import gzip
import os
import sys
from joblib import Parallel, delayed

def parse_sam(basename):
    with gzip.open(basename + '.sam.gz', 'rt') as sam:
        
        i = 0
        n = 0
        
        prot_counts = pd.Series()
        prot_counts.name = basename

        for line in sam:
            
            if line.startswith('@') == False:
                i = i + 1
                
                line = line.split('\t')
                    
                ## For mapped reads, add to tally for the reference sequence.  These are selected based on
                ## bitwise flag in the SAM file format and should only include the primary alignment:
                ## 1: template having multiple segments in sequencing
                ## 2: each segment properly aligned according to the aligner
                ## 16: SEQ being reverse complemented
                
                if line[1] not in ['77', '141']:            
                    rname = line[2]
                    n = n + 1
                    
                    try:
                        prot_counts[rname] = prot_counts[rname] + 1
                    except KeyError:
                        prot_counts[rname] = 1
        
                print('tallying hits for', basename + ':', 'hits =', n, 'out of', i, 'flag =', line[1]) 
        
        prot_counts.to_csv(basename + '.protcounts.txt')
        
basenames = []
        
for f in os.listdir('.'):
    if f.endswith('sam.gz'):            
            basename = f.split('.sam.gz')[0] 
            basenames.append(basename)
            
Parallel(n_jobs = -1, verbose = 5)(delayed(parse_sam)
(basename) for basename in basenames)

allprotcounts = pd.DataFrame()

for b in basenames:
    temp  = pd.read_csv(b + '.protcounts.txt', index_col = 0)
    allprotcounts = pd.concat([allprotcounts, temp], axis = 1)

allprotcounts.to_csv('allprotcounts.csv')