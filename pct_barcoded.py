#!/usr/bin/env python

'''calculate the percent of reads in a lane properly resolved by barcode
'''

import os,sys
from radtag_denovo import preprocess_radtag_lane
from glob import glob

fastq,analysis_folder = sys.argv[1:]

tot_reads = preprocess_radtag_lane.get_read_count(fastq)

indiv_barcoded = {}
fqs = glob('%s/*1_sequence.33.fq4*' % analysis_folder)
for fq in fqs:
    indiv_barcoded[fq] = preprocess_radtag_lane.get_read_count(fq)

print float(sum(indiv_barcoded.values()))/tot_reads
