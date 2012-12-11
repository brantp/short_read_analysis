#!/usr/bin/env python

'''
given a reference index file (.fai), a read group, a sample label and a .bam,

produces a sorted, RG tagged .bam

add_rg_and_sort_bam_by_ref.py ref.fai bam rg sm 
'''

import os, sys
from subprocess import Popen, PIPE

ref,bf,id,sm = sys.argv[1:]

headline = '@RG\tID:%s\tPL:Illumina\tLB:%s\tSM:%s\n' % (id,sm,sm)
ofh = open(bf[:-3]+'rg.sam','w')
print >> sys.stderr, 'adding RG to %s, output in %s' % (bf,bf[:-3]+'rg.sam')
ofh.write(headline)
for l in Popen('samtools view %s' % bf, shell=True,stdout=PIPE).stdout:
    ofh.write(l.strip('\n')+'\tRG:Z:'+id+'\n')
ofh.close()

os.system('samtools import %s %s %s' % (ref,bf[:-3]+'rg.sam',bf[:-3]+'rg_refhead.bam'))
os.system('samtools sort %s %s' % (bf[:-3]+'rg_refhead.bam', bf[:-3]+'rg_refsort'))
os.system('samtools index %s' % bf[:-3]+'rg_refsort.bam')

os.unlink(bf[:-3]+'rg.sam')
os.unlink(bf[:-3]+'rg_refhead.bam')
