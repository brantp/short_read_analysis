#!/usr/bin/env python

queue = 'general'
max_job_duration = 1440
job_ram = (2+1)*1024
job_batches = 500

import os,sys,SLURM,run_safe

geno,pheno,runs = sys.argv[1:]

basedir,basename = os.path.split(pheno)
donedir = os.path.join(basedir,os.path.splitext(basename)[0]+'-permute-donedir/')
logfile = os.path.join(basedir,os.path.splitext(basename)[0]+'-permute-logs/log-')

if not os.path.exists(donedir): os.makedirs(donedir)

trd = {}
for i in range(int(runs)):
    run_safe.add_cmd(trd,donedir+str(i),'random_wigs.py %s %s %s' % (geno,pheno,i) ,force_write=True)


#LSF.lsf_run_until_done(trd,logfile,queue,'','random-wigs',1000,3)
SLURM.run_until_done(trd,'random-wigs',logfile,max_job_duration,job_ram,job_batches,queue,MAX_RETRY=3)
