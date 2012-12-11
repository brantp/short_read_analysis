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
import os, sys, LSF
from glob import glob
from subprocess import Popen, PIPE
from collections import defaultdict
from short_read_analysis import preprocess_radtag_lane

RAM = 8 # gb ram for GATK
gatk_jar = '/n/hoekstra_lab/sureselect/Sting/dist/GenomeAnalysisTK.jar'
njobs = 40 

if __name__ == '__main__':

    outvcf, bwa_argstr, reference_fasta = sys.argv[1:4]
    reads = sys.argv[4:]

    outbase = os.path.dirname(outvcf)
    
    t = reference_fasta

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
    
    readbases = {}.fromkeys(list(set([rb[:-1 * (len('_sequence.33.fq4') + 1)] for rb in reads if rb.endswith('_sequence.33.fq4')])),'')
    readbases.update({}.fromkeys(list(set([rb[:-1 * (len('_sequence.33.fq4.gz') + 1)] for rb in reads if rb.endswith('_sequence.33.fq4.gz')])),'.gz'))

    #need to add support for _indexN.33.fq4 here
    
    for readbase,readext in readbases.items():

        readroot,rb = os.path.split(readbase)

        fc = os.path.basename(readroot)
        if rb.startswith('pass'): continue
        
        samf = os.path.join(readroot,'%spair-%s_bwa%s.sam' % (rb,tb,bp))
        sortbam = samf[:-3] + 'sort.bam'
        rg_ref_bam = samf[:-3] + 'rg_refsort.bam'
        rg_ref_bams.append(rg_ref_bam)
        #flowcell_by_bamstr[sortbam] = os.path.dirname(readroot)
        

        #bams.append(sortbam)
        ID = rb+fc
        SM = rb.split('_')[0]

        headline = '@RG\tID:%s\tPL:Illumina\tLB:%s\tSM:%s\n' % (ID,SM,SM)
        headfile = rg_ref_bam+'.headline'
        open(headfile,'w').write(headline)


        if os.path.exists(rg_ref_bam): continue #skip this individual if already mapped under these params

        r1 = os.path.join(readroot,'%s1_sequence.33.fq4%s' % (rb,readext))
        r2 = os.path.join(readroot,'%s2_sequence.33.fq4%s' % (rb,readext))
        sai1 = os.path.join(readroot,'%s1_sequence.33-%s_bwa%s.sai%s' % (rb,tb,bp,readext))
        sai2 = os.path.join(readroot,'%s2_sequence.33-%s_bwa%s.sai%s' % (rb,tb,bp,readext))
        print >> sys.stderr, 'aligning %s and %s:\n\tsam as %s, bam as %s' % (r1,r2,samf,sortbam)
        cmdstr = 'bwa aln %s %s %s -f %s; ' % (bwa_argstr,t,r1,sai1)
        cmdstr += 'bwa aln %s %s %s -f %s; ' % (bwa_argstr,t,r2,sai2)
        cmdstr += 'bwa sampe %s %s %s %s %s | ' % (t,sai1,sai2,r1,r2)
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


    if len(cmds) > 0:
        logfile = os.path.join(outbase,'bwa-%s-%s-log' % (bp,tb))
        jobids,namedict = LSF.lsf_jobs_submit(cmds,logfile,'normal_serial',jobname_base='bwa',num_batches=njobs)
        LSF.lsf_wait_for_jobs(jobids,logfile,namedict=namedict)

        cmds = LSF.lsf_no_success_from_log(logfile)
        while len(cmds) > 0:
            jobids,namedict = LSF.lsf_jobs_submit(cmds,logfile,'normal_serial',jobname_base='bwa')
            LSF.lsf_wait_for_jobs(jobids,logfile,namedict=namedict)
            cmds = LSF.lsf_no_success_from_log(logfile)

    cmds = []

    #for bam in bams:
    #    fc = flowcell_by_bamstr[bam]
    #    bamroot = os.path.dirname(bam)
    #    rg_ref_bam = bam[:-3] + 'sort.rg_refsort.bam'
    #    if os.path.exists(rg_ref_bam): continue
    #    ID = bam.split('_pair')[0]+'_'+fc
    #    SM = ID.split('_')[0]
    #    cmds.append('add_rg_and_sort_bam_by_ref.py %s.fai %s %s %s' % (reference_fasta,bam,ID,SM))
    
    #logfile = os.path.join(bamroot,'rg-sort-log')
    #jobids,namedict = LSF.lsf_jobs_submit(cmds,logfile,'normal_serial',jobname_base='rg-sort',num_batches=100)
    #LSF.lsf_wait_for_jobs(jobids,logfile,namedict=namedict)


    
    #COMPILE BAMS FROM READBASES
    istr = ' '.join(['-I '+fn for fn in rg_ref_bams])

    gatkcmd = 'java -Xmx%sg -jar  %s %s -R %s -T UnifiedGenotyper -o %s -all_bases -genotype' % (RAM,gatk_jar,istr,reference_fasta,outvcf)
    print >> sys.stderr, 'running GATK:\n%s' % gatkcmd

    os.system(gatkcmd)
