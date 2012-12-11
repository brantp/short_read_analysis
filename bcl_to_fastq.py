#!/usr/bin/env python

import os,sys
from glob import glob

def mkdir(outroot):
    try:
        os.makedirs(outroot)
        print >> sys.stderr, 'output root %s created' % outroot
    except:
        print >> sys.stderr, 'output root %s already exists' % outroot


I2B_JAR = '/n/home08/brantp/src/illumina2bam/dist/illumina2bam.jar'
S2F_JAR = '/n/home08/brantp/src/picard_svn/trunk/dist/SamToFastq.jar'


intensities, lane, outroot = sys.argv[1:]

bamroot = os.path.join(outroot,'bam')
fastqroot = os.path.join(outroot,'fastq')

mkdir(outroot)
mkdir(bamroot)
mkdir(fastqroot)

# 1) convert bcl to bam using illumina2bam. 0.5 - 1hr
bamfile = '%s/lane%s.bam' % (bamroot, lane)
if os.path.exists(bamfile):
    print >> sys.stderr, 'bam %s already exists' % bamfile
else:
    os.system('java -jar %s I=%s L=%s O=%s PF=false' % (I2B_JAR, intensities, lane, bamfile))

# 2) swap sequence indices for "index1" "index2" etc. 1 - 2hr?

idxbamfile = '%s/lane%s_indexedRG.bam' % (bamroot, lane)
if os.path.exists(idxbamfile):
    print >> sys.stderr, 'bam %s already exists' % idxbamfile
else:
    os.system('bam_RG_from_index.py %s | samtools view -Sb -o %s -' % (bamfile, idxbamfile))

# 3) split into fastq using picard (OPRG makes one fastq per read group, which is one per index)
os.system('java -Xmx2g -jar %s I=%s/lane%s_indexedRG.bam OPRG=true NON_PF=true ODIR=%s' % (S2F_JAR, bamroot, lane, fastqroot))
os.system('gzip -9 %s/*.fastq' % fastqroot)

fastqs = glob(os.path.join(fastqroot,'index*.fastq.gz'))
for f in fastqs:
    root,base = os.path.split(f)
    idxstr,read = base.split('.')[0].split('_')
    newf = os.path.join(root, 's_%s_%s_sequence_%s.txt.gz' % (lane,read,idxstr))
    print >> sys.stderr, 'rename %s to %s' % (f,newf)
    os.system('mv %s %s' % (f,newf))


