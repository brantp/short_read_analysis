#!/usr/bin/env python
'''quickie script drops terminal "B" (2) quality bases off ends of reads, retaining only reads longer than specified cutoff

works on one-line fastq
'''

import os,sys

maxlines = None

minlen = int(sys.argv[1])
fastq = sys.argv[2:]

outfns = [os.path.splitext(f)[0]+'-trimBs.txt' for f in fastq]
outfhs = [open(f,'w') for f in outfns]
infhs = [open(f) for f in fastq]

i=0
while 1:
    i+=1
    if maxlines is not None and i > maxlines:
        break
    try:
        lines = [fh.next() for fh in infhs]
        fields = [l.strip().rsplit(':',2) for l in lines]
        nonBquals = [f[-1].rstrip('B') for f in fields]
        numgood = [len(q) for q in nonBquals]
        if all([n >= minlen for n in numgood]):
            for f,n,fh in zip(fields,numgood,outfhs):
                fh.write('%s:%s:%s\n' % (f[0],f[1][:n],f[2][:n]))
    except StopIteration:
        break
    
[fh.close() for fh in outfhs]