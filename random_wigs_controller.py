#!/usr/bin/env python

queue = 'normal_serial'

import os,sys,LSF,run_safe

geno,pheno,runs = sys.argv[1:]

basedir = os.path.dirname(pheno)
donedir = os.path.join(basedir,'donedir/')
logfile = = os.path.join(basedir,'logs/log-')

if not os.path.exists(donedir): os.makedirs(donedir)

trd = {}
for i in range(int(runs)):
    run_safe.add_cmd(trd,donedir+str(i),'random_wigs.py %s %s %s' % (geno,pheno,i) ,force_write=True)


LSF.lsf_run_until_done(trd,logfile,queue,'','random-wigs',1000,3)
