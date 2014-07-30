#!/usr/bin/env python

'''
samtools_cmd_by_sample.py <bamfile> <samtools cmd> <outdir> <compress>

runs specified samtools command on bam subset for each sample in bamfile, output in outdir

if compress is one of [ gzip, bgzip ] the indicated compression will be used on the output file

(uses get_bam_rg_by_sample.py for extracting sample subsets; py_util/SLURM.py for parallellism)
'''

MAX_DUR = 60 #1 hour
MAX_RAM = 2*1024 #2GB ram
PARTITION = 'general' #SLURM queue 'general'

import SLURM, run_safe
import os,sys

from glob import glob

bamfile, samtools_cmd, outroot, compress = sys.argv[1:]

if compress in ('gzip','bgzip'):
    compress_ext = '.%s.gz' % compress
    compress_cmd = ' | %s -c ' % compress
else:
    compress_ext = ''
    compress_cmd = ''
    

#First off, get read group files
rgcmd = 'get_bam_rg_by_sample.py %s %s' % (bamfile, outroot)
ret = os.system(rgcmd)

if ret != 0:
    print >> sys.stderr, 'failed:\n%s' % rgcmd
    raise OSError

sample_files = glob(os.path.join(outroot,'*.rgids.txt'))
bambase = os.path.basename(bamfile)
stcmdstr = samtools_cmd.replace(' ','_')

trd = {}

for sf in sample_files:
    sm = os.path.basename(sf).rsplit('.',2)[0]
    outfile = os.path.join(outroot, '%s-%s-%s%s' % (bambase, sm, stcmdstr, compress_ext))
    cmd = 'samtools view -hR %s %s | samtools view -bS - | samtools %s /dev/stdin %s> %s' % (sf, bamfile, samtools_cmd, compress_cmd, outfile)
    slurmbase = os.path.join(outroot,'%s-%s-%s' % (bambase, sm, stcmdstr))
    run_safe.add_cmd(trd, slurmbase, cmd,force_write=True)


logfile = os.path.join(outroot,'%s-%s_SLURMwd/log' % (bambase,stcmdstr))
SLURM.run_until_done(trd,'samtools_by_indiv',logfile,MAX_DUR,MAX_RAM,100,PARTITION,MAX_RETRY=3)
