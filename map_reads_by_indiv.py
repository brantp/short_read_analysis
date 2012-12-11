#!/usr/bin/env python
'''
usage:
map_sureselect_lanes.py output.vcf "bwa arguments and values" reference_fasta indiv1.33.fq4 indiv2.33.fq4 ... indivN.33.fq4

an example of the bwa argument string might be "-l 40 -k 3 -q 3 -R 20"

resulting files will have a summary of arguments attached, e.g. -
indiv1-reference_fasta_bwa-l40-k3-q3-R20.sort.bam

all files except .vcf will be created in the same directory as the fastq data files
.vcf will be created as specified (e.g. specify full path to vcf!)


'''
import os, sys, re, LSF
from glob import glob
from subprocess import Popen, PIPE
from collections import defaultdict
from short_read_analysis import preprocess_radtag_lane

RAM = 16 # gb ram for GATK
gatk_jar = '/n/home08/brantp/src/GATK-git/dist/GenomeAnalysisTK.jar'
njobs = 40 

def find_paired_reads(reads,verbose=True):
    '''given a list of read filenames

    returns two lists (unpaired,paired) where unpaired is all single read files from reads
    paired is a list of 2-tuples (read1,read2) where read1 and read2 are fwd and rev reads
    '''

    paired = []
    unpaired = []
    skip = []
    for read in reads:
        if read in skip:
            continue
        
        readpath,readbase = os.path.split(read)
        if '_s_' in readbase:
            ind,voidstr,lane,suf = readbase.split('_',3)
            voidstr = '_'+voidstr #hack to get _s_ into mate filename
        else:
            voidstr = ''
            ind,lane,suf = readbase.split('_',2)

        if suf[0] == '1':
            #check for mate as read 2
            matebase = '_'.join([ind+voidstr,lane,'2'+suf[1:]])
            mate = os.path.join(readpath,matebase)
            if mate in reads:
                paired.append((read,mate))
                skip.append(read)
                skip.append(mate)
            else:
                if verbose:
                    print >> sys.stderr, '%s indicates read 1, \n\tno read 2 (%s) found; treat as unpaired' % (read,mate)
                unpaired.append(read)
                skip.append(read)
        elif suf[0] == '2':
            #check for mate as read 1
            matebase = '_'.join([ind,lane,'1'+suf[1:]])
            mate = os.path.join(readpath,matebase)
            if mate in reads:
                paired.append((mate,read))
                skip.append(read)
                skip.append(mate)
            else:
                if verbose:
                    print >> sys.stderr, '%s indicates read 2, no read 1 found; skip (use --unpair_reads to force inclusion)' % read

        else:
            #unpaired
            unpaired.append(read)
            skip.append(read)

    return unpaired,paired

if __name__ == '__main__':

    import argparse

    ds =  ' [%(default)s]'
    #create command line parser
    parser = argparse.ArgumentParser(description='performs BWA mapping for SINGLE READ data (all data will be interpreted as single read)')
    parser.add_argument('-v','--vcfname',default=None,help='if vcfname (filename) is supplied, will run GATK with --gatk_argstr and any additional BAM files specified using -a. VCF will be created in OUTROOT'+ds)
    parser.add_argument('-b','--bwa_argstr',default="'-q 3'",type=eval,help='arguments passed to BWA. Must be single AND double quoted for spaces'+ds)
    parser.add_argument('-g','--gatk_argstr',default="'-all_bases -genotype'",type=eval,help='arguments passed to GATK (only relevant if --outvcf specified) Must be single AND double quoted for spaces.'+ds)
    parser.add_argument('-a','--add_bam',default=[],action='append',help='additional BAM files for GATK (only relevant if --outvcf specified) any number of -a arguments may be supplied'+ds)
    parser.add_argument('-uc','--uncompressed_outputs',action='store_true',help='BWA outputs will be gzipped if not set'+ds)
    parser.add_argument('-up','--unpair_reads',action='store_true',help='will NOT attempt to pair read1 and read2 files in input reads'+ds)
    parser.add_argument('--debug',action='store_true',help='print commands but do not run anything'+ds)

    parser.add_argument('reference_fasta',help='reference for BWA')
    parser.add_argument('outroot',help='directory for logfile and vcf creation')
    parser.add_argument('reads',nargs='+',help='fastq for BWA (qualities should be fastq-33)')
    
    opts = parser.parse_args()

    vcfname, bwa_argstr, reference_fasta = opts.vcfname,opts.bwa_argstr,opts.reference_fasta
    reads = opts.reads

    outroot = opts.outroot
    
    t = reference_fasta
    if not os.path.exists(outroot): os.makedirs(outroot)

    tb = os.path.splitext(os.path.basename(t))[0]
    bp = bwa_argstr.replace(' ','')

    fai = t+'.fai'
    bwt = t+'.bwt'

    if not os.path.exists(fai): os.system('samtools faidx '+t)
    if not os.path.exists(bwt): os.system('bwa index '+t)
    
    cmds = []
    bams = []
    rg_ref_bams = []
    flowcell_by_bamstr = {}


    
    #readbases = {}.fromkeys(list(set([rb[:-1 * (len('_sequence.33.fq4') - 1)] for rb in reads if rb.endswith('_sequence.33.fq4')])),'')
    #readbases.update({}.fromkeys(list(set([rb[:-1 * (len('_sequence.33.fq4.gz') - 1)] for rb in reads if rb.endswith('_sequence.33.fq4.gz')])),'.gz'))

    if opts.uncompressed_outputs:
        readext = ''
    else:
        readext = '.gz'
        
    if opts.unpair_reads:
        print >> sys.stderr, '%s fastq processed (no PE check)' % len(reads)
    else:
        print >> sys.stderr, 'PE check invoked on %s fastq' % len(reads)
        unpaired,paired = find_paired_reads(reads)
        print >> sys.stderr, 'found %s PE; %s SR' % (len(paired),len(unpaired))
        reads = unpaired+paired

    for readfile in reads:
        print >> sys.stderr, readfile
        #throughout, tuple indicates PE
        if isinstance(readfile,tuple):
            readroot,readfilebase = os.path.split(readfile[0])
        else:
            readroot,readfilebase = os.path.split(readfile)

        #rb = readfilebase.split('_',1)[0]
        try:
            rb = re.search('(.+?)_[s\d]_',readfilebase).groups()[0]
        except:
            print >> sys.stderr, readfilebase
            raise

        fc = os.path.basename(readroot)
        if rb.startswith('pass'): continue

        if isinstance(readfile,tuple):
            samf = os.path.join(readroot,'%s_pair-%s_bwa%s.sam' % (rb,tb,bp))
        else:
            samf = os.path.join(readroot,'%s-%s_bwa%s.sam' % (rb,tb,bp))

        sortbam = samf[:-3] + 'sort.bam'
        rg_ref_bam = samf[:-3] + 'rg_refsort.bam'
        rg_ref_bams.append(rg_ref_bam)
        #flowcell_by_bamstr[sortbam] = os.path.dirname(readroot)
        

        #bams.append(sortbam)
        ID = rb+'_'+fc
        SM = rb

        headline = '@RG\tID:%s\tPL:Illumina\tLB:%s\tSM:%s\n' % (ID,SM,SM)
        headfile = rg_ref_bam+'.headline'
        if not opts.debug:
            open(headfile,'w').write(headline)


        if os.path.exists(rg_ref_bam):
            print >> sys.stderr, 'skip',rg_ref_bam
            continue #skip this individual if already mapped under these params

        if isinstance(readfile,tuple):
            r1,r2 = readfile
            sai1 = os.path.join(readroot,'%s_1-%s_bwa%s.sai%s' % (rb,tb,bp,readext))
            sai2 = os.path.join(readroot,'%s_2-%s_bwa%s.sai%s' % (rb,tb,bp,readext))
            print >> sys.stderr, 'aligning %s and %s:\n\tsam as %s, bam as %s' % (r1,r2,samf,sortbam)
            cmdstr = 'bwa aln %s %s %s -f %s; ' % (bwa_argstr,t,r1,sai1)
            cmdstr += 'bwa aln %s %s %s -f %s; ' % (bwa_argstr,t,r2,sai2)
            cmdstr += 'bwa sampe %s %s %s %s %s | ' % (t,sai1,sai2,r1,r2)
            cmdstr += 'samtools view -S - | add_rg_to_sam.py %s %s | ' % (ID,SM)
            cmdstr += 'samtools view -Sh -t %s.fai - | cat %s - | ' % (t, headfile)
            cmdstr += 'samtools view -bSh - | samtools sort - %s; ' % (rg_ref_bam[:-4])
            cmdstr += 'samtools index %s' % (rg_ref_bam)
        else:
            r1 = readfile
            sai1 = os.path.join(readroot,'%s-%s_bwa%s.sai%s' % (rb,tb,bp,readext))
            print >> sys.stderr, 'aligning %s:\n\tsam as %s, bam as %s' % (r1,samf,sortbam)
            cmdstr = 'bwa aln %s %s %s -f %s; ' % (bwa_argstr,t,r1,sai1)
            cmdstr += 'bwa samse %s %s %s | ' % (t,sai1,r1)
            cmdstr += 'samtools view -S - | add_rg_to_sam.py %s %s | ' % (ID,SM)
            cmdstr += 'samtools view -Sh -t %s.fai - | cat %s - | ' % (t, headfile)
            cmdstr += 'samtools view -bSh - | samtools sort - %s; ' % (rg_ref_bam[:-4])
            cmdstr += 'samtools index %s' % (rg_ref_bam)

        ##### old
        #cmdstr += 'bwa sampe %s %s %s %s %s > %s; ' % (t,sai1,sai2,r1,r2,samf )
        #cmdstr += 'samtools view -bS %s | samtools sort - %ssort; ' % (samf,samf[:-3])

        #cmdstr += 'add_rg_and_sort_bam_by_ref.py %s.fai %s %s %s' % (reference_fasta,bam,ID,SM)
        #####

        print >> sys.stderr, 'add:\n',cmdstr
    
        cmds.append(cmdstr)

    if opts.debug:
        print 'COMMANDS FOLLOW:\n'+'\n'.join(cmds)
    else:
        if len(cmds) > 0:
            import time
            logfile = os.path.join(outroot,'bwa-%s-%s-log' % (bp,tb))
            jobids,namedict = LSF.lsf_jobs_submit(cmds,logfile,'normal_serial',jobname_base='bwa',num_batches=njobs)
            time.sleep(20)
            LSF.lsf_wait_for_jobs(jobids,logfile,namedict=namedict)

            cmds = LSF.lsf_no_success_from_log(logfile)
            while len(cmds) > 0:
                jobids,namedict = LSF.lsf_jobs_submit(cmds,logfile,'normal_serial',jobname_base='bwa')
                time.sleep(20)
                LSF.lsf_wait_for_jobs(jobids,logfile,namedict=namedict)
                cmds = LSF.lsf_no_success_from_log(logfile)

            cmds = []

    
    #COMPILE BAMS FROM READBASES
    if vcfname is not None:
        outvcf = os.path.join(outroot,vcfname)
        istr = ' '.join(['-I '+fn for fn in list(set(rg_ref_bams+opts.add_bam))])
        
        gatkcmd = 'java -Xmx%sg -jar  %s %s -R %s -T UnifiedGenotyper -o %s %s' % (RAM,gatk_jar,istr,reference_fasta,outvcf,opts.gatk_argstr)
        print >> sys.stderr, 'running GATK:\n%s' % gatkcmd

        if not opts.debug:
            os.system(gatkcmd)
