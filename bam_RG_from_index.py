#!/usr/bin/env python

import os,sys,re

bamfile = sys.argv[1]

from subprocess import Popen, PIPE
from short_read_analysis.merge_read_Sample_by_index import idx_lookup
from radtag_denovo.preprocess_radtag_lane import match_index

os.system('samtools view -H '+bamfile)
for k,v in idx_lookup.items():
    print '@RG\tID:index%s\tPL:Illumina\tLB:%s\tSM:%s' % (v,k,v)
print '@RG\tID:PASS\tPL:Illumina\tLB:unknown\tSM:unknown'

inlinefeed = Popen('samtools view '+bamfile,shell=True,stdout=PIPE)

for l in inlinefeed.stdout:
    rtmatch = re.search('RT:Z:([ACGTN]+)',l)
    if rtmatch is None:
        print >> sys.stderr, 'line failed RT:Z: tag (readgroup)\n%s' % l
        raise ValueError
    idx = match_index(rtmatch.groups()[0],idx_lookup)
    if idx is None:
        print re.sub('RG:Z:.+?\t','RG:Z:PASS\t',l),
    else:
        print re.sub('RG:Z:.+?\t','RG:Z:index%s\t' % idx,l),
        
    
#need:
#@RG	ID:Egl16_2_101104_B80EP1ABXX_lane2	PL:Illumina	LB:Egl16	SM:Egl16
