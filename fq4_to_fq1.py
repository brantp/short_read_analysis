#!/usr/bin/env python
'''
convert 4 line fastq to 1 line fastq
'''

import os, sys
infile = sys.argv[1]

qnext = False
snext = False
q = ''
s = ''
h = ''

for l in open(infile):
    if snext:
        s = l.strip()
        snext = False
    elif qnext:
        q = l
        print ':'.join([h,s,q]),
        qnext = False
    elif l.startswith('@'):
        h = l[1:].strip()
        snext = True
    elif l.startswith('+'):
        qnext = True
