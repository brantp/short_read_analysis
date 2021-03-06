#!/usr/bin/env python

'''
get_bam_rg_by_sample.py <bam> <outdir>

takes a bam file and extracts RG header lines, groups by sample and writes sample-labeled files to <outdir>
'''

import os, sys
from subprocess import Popen, PIPE
from collections import defaultdict

bamfile, outroot = sys.argv[1:]

try:
    os.makedirs(outroot)
except:
    pass

header_lines = Popen('samtools view -H %s | grep -P "^@RG" ' % bamfile, shell = True, stdout = PIPE).stdout.readlines()

id_by_sm = defaultdict(list)

for l in header_lines:
    rg_d = dict([el.split(':') for el in l.strip().split()[1:]])
    id_by_sm[ rg_d['SM'] ].append( rg_d['ID'] )

for k,v in id_by_sm.items():
    outf = os.path.join(outroot,k+'.rgids.txt')
    outfh = open(outf,'w')
    for rgid in v:
        outfh.write(rgid+'\n')
    outfh.close()
