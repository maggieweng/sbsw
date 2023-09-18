# -*- coding: utf-8 -*-
"""
Created on Sat Aug 14 07:47:03 2021

@author: jeff
"""

import pandas as pd
import gzip
import os
import sys

name = sys.argv[1]
cwd = os.getcwd() + '/' # The current working directory.

prot_counts = pd.Series()
prot_counts.name = 'n_hits'
    
i = 0
f = 0
    
with gzip.open(cwd + name + '.sam.gz', 'rt') as sam:
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
                f = f + 1
                
                try:
                    prot_counts[rname] = prot_counts[rname] + 1
                except KeyError:
                    prot_counts[rname] = 1
    
            print('tallying hits for', name + ':', 'hits =', f, 'out of', i, 'flag =', line[1])
            
## import GFF - only if using prokka annotation
            
#col_names = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']           
#annot_gff = pd.read_csv('KBase_derived_mygenome.gff_edited.gff', sep = '\t', index_col = 0, names = col_names)

## import custom pgap compilation - only if using PGAP
            
annot_gff = pd.read_csv('210924_bins.annotations.csv', index_col = 0)

## merge

prot_counts = prot_counts.reindex(annot_gff.index)
prot_unique_cds_df = pd.concat([annot_gff, prot_counts], axis = 1)
prot_unique_cds_df.dropna(subset = ['n_hits'], inplace = True)
prot_unique_cds_df.sort_values(by = 'n_hits', inplace = True, ascending = False)

prot_unique_cds_df.to_csv(name + '_tally.csv')