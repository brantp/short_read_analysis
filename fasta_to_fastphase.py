#!/usr/bin/env python

import os,sys
import Aln


fasta = sys.argv[1]

fastphase = os.path.splitext(fasta)[0]+'.inp'

alnfa = Aln.AlnFasta(fasta)

ofh = open(fastphase,'w')

ofh.write('%s\n%s\n' % (len(alnfa),len(alnfa.values()[0])))

for ind, s in alnfa.items():
    ofh.write('# %s\n%s\n%s\n' % ((ind,)+tuple([si.upper().replace('N','?') for si in s.undegenerate(include_n=False,return_both=True)])))

ofh.close()
