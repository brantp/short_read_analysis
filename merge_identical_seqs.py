#!/usr/bin/env python
'''read through fastq (1 line format) and return unique sequences in fasta

takes n_bases_to_skip fastq_file

output:
>num_found.median_qual
sequence

'''

import os,sys,re,Util,numpy

from collections import defaultdict

n_skip,fastq = sys.argv[1:]
n_skip = int(n_skip)

seqs = defaultdict(int)
quals = defaultdict(list)

for l in open(fastq):
    s,q = l.split(':')[-2:]
    s,q = s[n_skip:],[ord(c)-64 for c in q[n_skip:-1]]
    seqs[s] += 1
    quals[s].append(q)

for k,v in quals.items():
    arv = numpy.array(v)
    medv = [numpy.median(pos) for pos in arv.transpose()]
    medqstr = ''.join([chr(int(i)+64) for i in medv])
    print '>%s.%s\n%s' % (seqs[k],medqstr,k)
    
