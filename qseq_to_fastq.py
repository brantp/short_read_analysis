#!/usr/bin/env python
#
#convert qseq to one-line fastq

import sys

for l in open(sys.argv[1]):
    elem = l.split()
    elem[1] = int(elem[1])
    elem[8] = elem[8].replace('.','N')
    print '%s_%04d:%s:%s:%s:%s#%s/%s:%s:%s' % tuple(elem[:-1])
