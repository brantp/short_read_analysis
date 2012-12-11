#!/usr/bin/env python
'''given paired end fq1 or fq4 files
usage:
process_sureselect_lanes.py outroot s_lane_1_sequence.fq s_lane_2_sequence.fq

generates sanger 33-based 4-line fastq for each individual in outroot
'''

import os, sys, LSF, re
from glob import glob
from subprocess import Popen, PIPE
from collections import defaultdict
from radtag_denovo import preprocess_radtag_lane
from Util import smartopen as open

def get_read_count(filename,lnum):
    
    if filename.endswith('.gz'):
        print >> sys.stderr, 'getting read count for compressed file',filename,'...',
        rc = int(Popen('gunzip -c %s | wc -l' % filename,shell=True,stdout=PIPE).stdout.read().split()[0]) / lnum
        print >> sys.stderr, rc
        return rc
    else:
        print >> sys.stderr, 'getting read count for file',filename,'...',
        rc = int(Popen('wc -l %s' % filename,shell=True,stdout=PIPE).stdout.read().split()[0]) / lnum
        print >> sys.stderr, rc
        return rc

def generate_unfinished_seqret(outroot,rext,lnum):
    cmds = []
    for f in glob(os.path.join(outroot,'*'+rext)):
        if os.path.basename(f).startswith('pass'): continue
        print >> sys.stderr, 'consider',f,'...',
        qfh = open(f)
        fbaseQ = None
        while fbaseQ is None:
            fbaseQ = preprocess_radtag_lane.get_baseQ(preprocess_radtag_lane.next_read_from_fh(qfh,lnum)[2])
        qfh.close()
        if fbaseQ == 33:
            print >> sys.stderr, 'qualities already sanger-33 encoded, create .33.fq4 symlink.'
            if f.endswith('.gz'):
                f4 = f[:-6]+'33.fq4.gz'
            else:
                f4 = f[:-3]+'33.fq4'
            os.system('ln -s %s %s' % (f,f4))
            continue
        else:
            print >> sys.stderr, 'qualities illumina-64 encoded, queue conversion.'

        if f.endswith('.gz'):
            f4 = f[:-6]+'33.fq4.gz'
            if (not os.path.exists(f4)) or (os.path.getsize(f4) < 1000 and get_read_count(f4,lnum) == 0):
                cmds.append("gunzip -c %s | seqret fastq-illumina::stdin fastq::stdout | gzip > %s" % (f,f4))
                print >> sys.stderr, 'run compressed'
            else:
                print >> sys.stderr, 'already present'
        else:
            f4 = f[:-3]+'33.fq4'
            if (not os.path.exists(f4)) or (os.path.getsize(f4) < 1000 and get_read_count(f4,lnum) == 0):
                cmds.append("seqret fastq-illumina::%s fastq::%s" % (f,f4))
                print >> sys.stderr, 'run'
            else:
                print >> sys.stderr, 'already present'
                
                
    return cmds

if __name__ == '__main__':

    try:
        outroot,fread,rread = sys.argv[1:]
        paired=True
    except ValueError:
        outroot,fread = sys.argv[1:]
        paired=False

    if fread.endswith('.gz'):
        rext = os.path.splitext(fread[:-3])[-1] + '.gz'
    else:
        rext = os.path.splitext(fread)[-1]

    fh = open(fread)
    chr1 = fh.read(1)
    fh.close()

    if chr1 == '@':
        lnum = 4
        print >> sys.stderr, 'using 4-line fastq'
    else:
        lnum = 1
        print >> sys.stderr, 'using 1-line fastq'

    qfh = open(fread)
    baseQ = None
    while baseQ is None:
        baseQ = preprocess_radtag_lane.get_baseQ(preprocess_radtag_lane.next_read_from_fh(qfh,lnum)[2])
    qfh.close()

    split_done_file = os.path.join(outroot,'split_done')
    if not os.path.exists(split_done_file):

        nreads_fwd = get_read_count(fread,lnum)
        if paired:
            nreads_rev = get_read_count(rread,lnum)
            if nreads_fwd != nreads_rev:
                raise ValueError, 'read 1 and read 2 files differ in size; very bad things.'

        nticks = 20

        print >> sys.stderr, 'splitting reads for individuals'
        indiv_data = preprocess_radtag_lane.get_individual_data_for_lane(fread)
        nreads = nreads_fwd
        if paired:
            fqs = (fread,rread)
            fqfhs = [open(f) for f in fqs]
            irop = (os.path.join(outroot,'%'+re.sub('index\d+','sequence',os.path.basename(fread))),os.path.join(outroot,'%'+re.sub('index\d+','sequence',os.path.basename(rread))))
            passfh = [open(p % 'pass','w') for p in irop]
        else:
            fqs = fread
            fqfhs = open(fqs)
            irop = os.path.join(outroot,'%'+re.sub('index\d+','sequence',os.path.basename(fread)))
            passfh = open(irop % 'pass','w')
        
        
        fhdict = {}
        i = 0
        tickon = nreads/nticks
        print >> sys.stderr, '\tloading'

        for i in xrange(nreads):
            if i%tickon==0: print >> sys.stderr, '\t\t%s / %s' % (i,nreads)
            if paired:
                lines = tuple([preprocess_radtag_lane.next_read_from_fh(fh,lnum) for fh in fqfhs])
                indiv,read,qual = preprocess_radtag_lane.assign_read_to_indiv(lines,indiv_data,indiv_reads_out_pattern=irop,fhdict=fhdict,passfh=passfh,read2_has_idx=True,lnum=4,baseQ_in=baseQ)
            else:
                line = preprocess_radtag_lane.next_read_from_fh(fqfhs,lnum)
                indiv,read,qual = preprocess_radtag_lane.assign_read_to_indiv(line,indiv_data,indiv_reads_out_pattern=irop,fhdict=fhdict,passfh=passfh,lnum=4,baseQ_in=baseQ)

        if paired:
            for fh in passfh:
                fh.close()
        else:
            passfh.close()


        for fh in fhdict.values():
            fh.close()

        print >> sys.stderr, 'read splitting done'
        os.system('touch '+split_done_file)
    
    print rext

    cmds = generate_unfinished_seqret(outroot,rext,lnum)
    import time
    while len(cmds) > 0:
        logfile = os.path.join(outroot,'seqret-log')
        jobids,namedict = LSF.lsf_jobs_submit(cmds,logfile,'short_serial',jobname_base='seqret',num_batches=10)
        time.sleep(20)
        LSF.lsf_wait_for_jobs(jobids,logfile,namedict=namedict)
        cmds = generate_unfinished_seqret(outroot,rext,lnum)

    '''
    for cmd in cmds:
        print >> sys.stderr, cmd
        ret = os.system(cmd)
        if ret != 0:
            sys.exit()
    '''
