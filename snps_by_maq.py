#!/usr/bin/env python

'''given:
 settings dict file
 reference
 one or more reads files

follows maq snp pipeline per:
maq map
maq assemble
maq cns2snp
maq.pl SNPfilter 

performing necessary steps to generate binaries and interconvert formats as indicated by file extensions
'''

import os, sys, re

settings,reference = sys.argv[1:3]
reads = sys.argv[3:]
for rf in reads:
    base,ext = os.path.splitext(rf)
    print base,ext
    if ext == '.txt': #check for 4-line fastq, i.e. starts with '@'
    	firstchar = open(rf).read(1)
    	if firstchar == '@':
    	    print >>sys.stderr, 'converting %s to 4-line fastq %s' % #totally unfinished; probably dropping maq!