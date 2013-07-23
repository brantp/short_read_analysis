#!/usr/bin/env python

queue = 'normal_serial'

import os,sys,LSF,run_safe
from glob import glob

genodir,pheno = sys.argv[1:]

basedir = os.path.dirname(pheno)
#donedir = os.path.join(basedir,'donedir/')
logfile = os.path.join(basedir,'logs/log-')

#if not os.path.exists(donedir): os.makedirs(donedir)

genos = glob(os.path.join(genodir,'*-geno.txt'))

print >> sys.stderr, 'wigs on %s contigs' % len(genos)

trd = {}
for geno in genos:
    outbase = geno.replace('-geno.txt','_output-')
    run_safe.add_cmd(trd,outbase[:-1],'wigs simevo -g %s -s %s -f %s -x 40000 -b 500 -t 100 -d 1' % (geno,pheno,outbase) ,force_write=True)

LSF.lsf_run_until_done(trd,logfile,queue,'','wigs-by-chrom',500,3)
