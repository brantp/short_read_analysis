#!/usr/bin/env python
'''
usage:
map_reads_by_indiv-stampy.py output.vcf "stampy arguments and values" reference_fasta indiv1.33.fq4 indiv2.33.fq4 ... indivN.33.fq4

an example of the stampy argument string might be "--substitutionrate=0.05"

resulting files will have a summary of arguments attached, e.g. -
indiv1-reference_fasta_stampy--substitutionrate_0.05.sort.bam

all files except .vcf will be created in the same directory as the fastq data files
.vcf and LSF logs will be created in <outroot> 


'''
import os, sys, re, LSF, numpy
from glob import glob
from subprocess import Popen, PIPE
from collections import defaultdict
from radtag_denovo import preprocess_radtag_lane
import run_safe
import time

compress_vcf = True
picardRAM = 4
#gatk_jar = '/n/home08/brantp/src/GATK-git/dist/GenomeAnalysisTK.jar'
#gatk2_jar = '/n/home08/brantp/src/GenomeAnalysisTK-2.1-8-g5efb575/GenomeAnalysisTK.jar'
#gatk2_jar = '/n/home08/brantp/src/GenomeAnalysisTK-2.2-16-g9f648cb/GenomeAnalysisTK.jar'
gatk2_jar = '/n/home08/brantp/src/GenomeAnalysisTK-2.4-3-g2a7af43/GenomeAnalysisTK.jar'
gatk_jar = gatk2_jar
#picard_jar_root = '/n/home08/brantp/src/picard_svn/trunk/dist'
picard_jar_root = '/n/home08/brantp/src/picard_svn_20130220/trunk/dist'
picard_jar = os.path.join(picard_jar_root,'MergeSamFiles.jar')
picard_seqdict_jar = os.path.join(picard_jar_root,'CreateSequenceDictionary.jar')
#stampy_module = 'bio/stampy-1.0.18'
stampy_module = 'bio/stampy-1.0.21_python2.7.3'
min_ind_realign = 12
MAX_RETRY = 3
MERGE_BAMS_ABOVE = 50

skip_contigs = [] #['ruf_bac_7180000001736']

def unfinished_cmds(to_run_dict,finished_ext='.done'):
    cmds = []
    for finished_base,cmd in to_run_dict.items():
        if not os.path.exists(finished_base+finished_ext):
            cmds.append(cmd)
    return cmds

def idxstr_from_idx(x):
    return x and '_index%s' % x or ''

def seq_len_from_fasta(fasta,return_order=False):
    dictfile = os.path.splitext(fasta)[0]+'.dict'
    if not os.path.exists(dictfile):
        os.system('java -jar %s REFERENCE=%s OUTPUT=%s' % (picard_seqdict_jar,fasta,dictfile))
    seq_len = {}

    if return_order:
        order = []
        for l in open(fasta):
            if l.startswith('>'):
                order.append(l[1:].strip().split()[0])

    for l in open(dictfile):
        if l.startswith('@SQ'):
            f = l.strip().split()
            seq_len[f[1][3:]] = int(f[2][3:])

    if return_order:
        return seq_len,order
    else:
        return seq_len

def partition_reference(fasta,parts,include_regions=None):

    seq_len,order = seq_len_from_fasta(fasta,True)
    if include_regions is None:
        partlen = numpy.ceil(sum(seq_len.values())/float(parts))
    else:
        partlen = numpy.ceil(sum([v for k,v in seq_len.items() if k in include_regions])/float(parts))
    run_parts = [[]]
    incl = 0
    for s in order:
        l = seq_len[s]
        if include_regions is not None and not s in include_regions:
            continue
        seqpos = 0
        while seqpos < l:
            if (l-seqpos) + incl < partlen:
                run_parts[-1].append('%s:%s-%s' % (s,1+int(seqpos),int(l)))
                incl += (l-seqpos)
                seqpos += (l-seqpos)
            else:
                incl_this = partlen-incl
                incl_to = seqpos+incl_this
                run_parts[-1].append('%s:%s-%s' % (s,1+int(seqpos),int(incl_to)))
                run_parts.append([])
                incl = 0
                seqpos += incl_this
    return run_parts
        
def realign_bams_lsf(bams,ref,outroot,njobs,min_ind_realign,queue='normal_serial',job_ram='20000',targetcreator_opts='--maxIntervalSize 5000',gatk_jar=gatk_jar,gatk_ram=8,force_links=False,MAX_RETRY=MAX_RETRY,fallback_queue=''):
    '''force_links replaces existing symlinks
    job ram in MB
    '''
    realign_root = os.path.join(outroot,'realign')
    intervals_parts_root = os.path.join(realign_root,'parts')
    intervals_file = os.path.join(intervals_parts_root,'all_part.RealignTargets.intervals')
    try:
        os.makedirs(intervals_parts_root)
    except:
        pass


    os.chdir(realign_root)
    bams_to_link = {}
    link_to_realign = {}
    realign_to_srcroot = {}
    bams_to_workdir = {}
    realigned_bams = []
    
    for bam in bams:
        bamsrcroot,bambase = os.path.split(bam)
        fcroot = os.path.basename(bamsrcroot)
        bamlink = os.path.join(realign_root,fcroot,bambase)
        realigned = bamlink[:-4]+'.realigned.bam'
        realignedlink = bam[:-4]+'.realigned.bam'
        realignedidx = bamlink[:-4]+'.realigned.bai'
        realignedlinkidx = bam[:-4]+'.realigned.bai'
        
        realigned_bams.append(realignedlink)
        if force_links:
            try:
                os.unlink(realignedlink)
            except:
                pass
            try:
                os.unlink(realignedlinkidx)
            except:
                pass
        if not os.path.exists(realignedlink):
            bams_to_link[bam] = bamlink
            link_to_realign[bamlink] = realigned
            realign_to_srcroot[realigned] = realignedlink
            bams_to_workdir[bam] = os.path.dirname(bamlink)
            if not os.path.exists(bamlink):
                try:
                    os.makedirs(os.path.dirname(bamlink))
                except:
                    pass
                ret = os.system('ln -s %s %s' % (bam,bamlink))
                if ret != 0:
                    raise OSError, 'ln failed: %s -> %s' % (bam,bamlink)

    if len(bams_to_link) == 0: #nothing to do here; return links
        print >> sys.stderr, 'realigned bams all present'
        return realigned_bams

    #otherwise get down to business
    print >> sys.stderr, 'Perform realignment:'
    if os.path.exists(intervals_file+'.done'):
        print >> sys.stderr, 'using %s as intervals_file' % intervals_file
    else:
        bamstr = ' -I '.join(bams)
        intervals_parts_regions = partition_reference(ref,njobs)
        to_run_dict = {}
        for i,part in enumerate(intervals_parts_regions):
            reg_str = ' -L '.join(part)
            intervals_parts_file = os.path.join(intervals_parts_root,'part%s.RealignTargets.intervals' % (i))
            cmd = 'java -Xmx%sg -jar %s -T RealignerTargetCreator -I %s -R %s -L %s %s -o %s' % (gatk_ram,gatk_jar,bamstr,ref,reg_str,targetcreator_opts,intervals_parts_file)
            to_run_dict[intervals_parts_file] = run_safe.safe_script(cmd,intervals_parts_file)

        logfile = os.path.join(intervals_parts_root,'logs','RealignerTargetCreator')

        #replace with run_until_done
        LSF.lsf_run_until_done(to_run_dict,logfile,queue,'-R "select[mem>%s]"' % job_ram, 'targetcreator',njobs,MAX_RETRY)
        if fallback_queue:
            LSF.lsf_run_until_done(to_run_dict,logfile,fallback_queue,'-R "select[mem>%s]"' % job_ram, 'targetcreator',njobs,MAX_RETRY)
        
        #cmds = unfinished_cmds(to_run_dict)
        #while cmds:
        #    jobids,namedict = LSF.lsf_jobs_submit(cmds,logfile,queue,'-R "select[mem>%s]"' % job_ram, jobname_base='targetcreator')
        #    LSF.lsf_wait_for_jobs(jobids,logfile,queue,namedict=namedict)
        #    cmds = unfinished_cmds(to_run_dict)

        catcmd = 'cat %s > %s' % (' '.join(to_run_dict.keys()),intervals_file)
        ret = os.system(run_safe.safe_script(catcmd,intervals_file))
        if ret != 0:
            raise OSError, 'cat failed: %s' % catcmd

            
    samples_per_batch = max(min_ind_realign, len(bams_to_link)/njobs)
    bam_batches_by_dir = []
    this_batch = []
    lastroot = None
    for bam in sorted(bams_to_link.keys()):
        if lastroot != bams_to_workdir[bam] or len(this_batch) == samples_per_batch:
            if this_batch:
                bam_batches_by_dir.append((lastroot,this_batch))
            this_batch = []
            lastroot = bams_to_workdir[bam]
        this_batch.append(bam)
    if this_batch: #process last batch
        bam_batches_by_dir.append((lastroot,this_batch))

    if len(bam_batches_by_dir) == 0:
        print >> sys.stderr, 'realignments present'
    else:
        print >> sys.stderr, 'REALIGNMENT BATCH SUMMARY:'
        for workdir,bam_batch in bam_batches_by_dir:
            print >> sys.stderr, '\tbams: %s working: %s' % (len(bam_batch),workdir)
            #for bam in bam_batch:
            #    print >> sys.stderr, '\t\t%s' % bam

        to_run_dict = {}
        for i,(workdir,bam_batch) in enumerate(bam_batches_by_dir):
            bamstr = ' -I '.join(bam_batch)
            donefile = os.path.join(realign_root,'realign_batch%sof%s' % (i,len(bam_batches_by_dir)))
            cmd = 'cd %s; java -Xmx%sg -jar %s -T IndelRealigner -model USE_SW -I %s -R %s --targetIntervals %s -nWayOut .realigned.bam' % (workdir,gatk_ram,gatk_jar,bamstr,ref,intervals_file)
            to_run_dict[donefile] = run_safe.safe_script(cmd,donefile)

        logfile = os.path.join(realign_root,'logs','IndelRealigner')
        LSF.lsf_run_until_done(to_run_dict,logfile,queue,'-R "select[mem>%s]"' % job_ram, 'realigner',njobs,MAX_RETRY)
        if fallback_queue:
            LSF.lsf_run_until_done(to_run_dict,logfile,fallback_queue,'-R "select[mem>%s]"' % job_ram, 'realigner',njobs,MAX_RETRY)
        #cmds = unfinished_cmds(to_run_dict)
        #while cmds:
        #    jobids,namedict = LSF.lsf_jobs_submit(cmds,logfile,queue,'-R "select[mem>%s]"' % job_ram, jobname_base='realigner')
        #    LSF.lsf_wait_for_jobs(jobids,logfile,queue,namedict=namedict)
        #    cmds = unfinished_cmds(to_run_dict)

    for realigned,realignedlink in realign_to_srcroot.items():
        if not os.path.exists(realignedlink):
            ret = os.system('ln -s %s %s' % (realigned,realignedlink))
            if ret != 0:
                raise OSError, 'ln failed: %s -> %s' % (realigned,realignedlink)
        if not os.path.exists(realignedlink[:-1]+'i'):
            ret = os.system('ln -s %s %s' % (realigned[:-1]+'i',realignedlink[:-1]+'i'))
            if ret != 0:
                raise OSError, 'ln failed: %s -> %s' % (realigned[:-1]+'i',realignedlink[:-1]+'i')

    return realigned_bams

def start_end_strs(li):
    c,b,e = re.search('(.+?):(\d+)-(\d+)',li[0]).groups()
    start = '%s-%s' % (c,b)
    c,b,e = re.search('(.+?):(\d+)-(\d+)',li[-1]).groups()
    end = '%s-%s' % (c,e)
    return start,end

def call_variants_gatk_lsf(bams,ref,outroot,vcfbase,njobs=100,gatk_program='UnifiedGenotyper',gatk_args='-out_mode EMIT_ALL_CONFIDENT_SITES -dcov 50 -glm BOTH',gatk_jar=gatk_jar,gatk_ram=4,tmpdir=None,queue='normal_serial',job_ram='30000',MAX_RETRY=MAX_RETRY,include_regions=None,compress_vcf=True,fallback_queue=''):
    if tmpdir is None:
        tmpdir = os.path.join(outroot,'gatk_tmp')
    bamstr = ' -I '.join(bams)
    regions = partition_reference(ref,njobs,include_regions)
    vcfbasename = vcfbase.endswith('.vcf') and vcfbase[:-4] or vcfbase
    gatkoutvcfbase = '%s-GATK-%s' % (vcfbasename,gatk_program)
    if compress_vcf:
        vcfext = '.vcf.gz'
    else:
        vcfext = '.vcf'
        
    gatkoutvcf = os.path.join(outroot,gatkoutvcfbase+vcfext)
    vcf_parts_root = os.path.join(outroot,gatkoutvcfbase+'-vcf_parts')
    try:
        os.makedirs(vcf_parts_root)
    except:
        pass

    to_run_dict = {}
    for i,reg in enumerate(regions):
        reg = [r for r in reg if not r.split(':')[0] in skip_contigs]
        if len(reg) == 0:
            continue
        start,end = start_end_strs(reg)
        regstr = ' -L '.join(reg)
        partvcf = os.path.join(vcf_parts_root,'%s_%dof%d_%sto%s%s' % (gatkoutvcfbase,i,len(regions),start,end,vcfext))
        part_sh = os.path.join(vcf_parts_root,'%s_%dof%d_%sto%s.sh' % (gatkoutvcfbase,i,len(regions),start,end))
        cmd = 'java -Xmx%sg -Djava.io.tmpdir=%s -jar  %s -R %s -T %s -o %s %s -I %s -L %s' % (gatk_ram,tmpdir,gatk_jar,ref,gatk_program,partvcf,gatk_args,bamstr,regstr)
        #open(part_sh,'w').write('#!/usr/bin/env bash\n'+cmd+'\n')
        #os.system('chmod +x %s' % part_sh)
        to_run_dict[partvcf] = run_safe.safe_script(cmd,partvcf,force_write=True)

    logfile = os.path.join(vcf_parts_root,'logs',gatk_program)
    LSF.lsf_run_until_done(to_run_dict,logfile,queue,'-R "select[mem>%s]"' % job_ram, 'gatk',njobs,MAX_RETRY)
    if fallback_queue:
        LSF.lsf_run_until_done(to_run_dict,logfile,fallback_queue,'-R "select[mem>%s]"' % job_ram, 'gatk',njobs,MAX_RETRY)
    #cmds = unfinished_cmds(to_run_dict)
    #while cmds:
    #    jobids,namedict = LSF.lsf_jobs_submit(cmds,logfile,'normal_serial','-R "select[mem>30000]"', jobname_base='gatk')
    #    LSF.lsf_wait_for_jobs(jobids,logfile,'normal_serial',namedict=namedict)
    #    cmds = unfinished_cmds(to_run_dict)

    #varstr = ' -V '.join(sorted(to_run_dict.keys()))
    #cmd = 'java -Xmx%sg -jar %s -T CombineVariants --assumeIdenticalSamples -R %s -V %s -o %s' % (gatk_ram,gatk_jar,ref,varstr,gatkoutvcf)
    cmd = merge_vcf_parts_cmd(to_run_dict.keys(),ref,gatkoutvcf,gatk_jar,gatk_ram,tmpdir)
    ret = os.system(run_safe.safe_script(cmd,gatkoutvcf))
    if ret != 0:
        raise OSError, 'VCF merge failed:\n%s' % cmd

def merge_vcf_parts_cmd(vcfparts,ref,outvcf,gatk_jar,gatk_ram,tmpdir,rod_type = ''):
    #varstr = (' -V%s ' % rod_type).join(sorted(vcfparts))
    #cmd = 'java -Xmx%sg -Djava.io.tmpdir=%s -jar %s -T CombineVariants --assumeIdenticalSamples -R %s -V %s -o %s' % (gatk_ram,tmpdir,gatk_jar,ref,varstr,outvcf)
    varstr = ' '.join(sorted(vcfparts))
    cmd = 'merge_vcf_parts.py -r %s -v %s -o %s' % (ref,varstr,outvcf)
    return cmd

def call_variants_mpileup_lsf(bams,ref,outroot,vcfbase,njobs=100,mpileup_args='',gatk_jar=gatk_jar,gatk_ram=8,tmpdir=None,queue='normal_serial',job_ram='30000',MAX_RETRY=MAX_RETRY,include_regions=None):
    if tmpdir is None:
        tmpdir = os.path.join(outroot,'gatk_tmp')

    bamstr = ' -I '.join(bams)
    regions = partition_reference(ref,njobs,include_regions)
    vcfbasename = vcfbase.endswith('.vcf') and vcfbase[:-4] or vcfbase
    mpoutvcfbase = '%s-mpileup' % (vcfbasename)
    mpoutvcf = os.path.join(outroot,mpoutvcfbase+'.vcf')
    vcf_parts_root = os.path.join(outroot,mpoutvcfbase+'-vcf_parts')
    try:
        os.makedirs(vcf_parts_root)
    except:
        pass

    to_run_dict = {}
    #merge_subparts_trd = {}
    subparts = []
    for i,reg in enumerate(regions):
        start,end = start_end_strs(reg)
        #regstr = ' -L '.join(reg)
        partdonebase = os.path.join(vcf_parts_root,'%s_%dof%d_%sto%s-parts' % (mpoutvcfbase,i,len(regions),start,end))
        partvcf = os.path.join(vcf_parts_root,'%s_%dof%d_%sto%s.vcf' % (mpoutvcfbase,i,len(regions),start,end))
        part_sh = os.path.join(vcf_parts_root,'%s_%dof%d_%sto%s.sh' % (mpoutvcfbase,i,len(regions),start,end))
        #cmd = 'java -Xmx%sg -Djava.io.tmpdir=%s -jar  %s -R %s -T %s -o %s %s -I %s -L %s' % (gatk_ram,tmpdir,gatk_jar,ref,gatk_program,partvcf,gatk_args,bamstr,regstr)

        this_trd = {}
        for this_reg in reg:
            subpart = this_reg.split(':')[0]
            subpartvcf = os.path.join(vcf_parts_root,'%s_%dof%d_%sto%s-%s.vcf' % (mpoutvcfbase,i,len(regions),start,end,subpart))
            this_cmd = 'samtools mpileup -Dgu -r %s %s -f %s %s | bcftools view -cvg  - > %s 2>  %s.log' % (this_reg,mpileup_args,ref,bamstr,subpartvcf,subpartvcf)
            this_trd[subpartvcf] = run_safe.safe_script(this_cmd,subpartvcf)
            
        cmd_parts = unfinished_cmds(this_trd)
        cmd = '; '.join(cmd_parts)
        #open(part_sh,'w').write('#!/usr/bin/env bash\n'+cmd+'\n')
        #os.system('chmod +x %s' % part_sh)
        to_run_dict[partdonebase] = run_safe.safe_script(cmd,partdonebase)
        subparts.extend(this_trd.keys())
        #vcfparts = this_trd.keys() ### <---MAKE THIS WORK (merge in parts before merge all)
        #merge_subparts_trd[partvcf] = run_safe.safe_script(merge_vcf_parts_cmd(vcfparts,ref,partvcf,gatk_jar,gatk_ram,tmpdir,rod_type = ':VCF'),partvcf)

    logfile = os.path.join(vcf_parts_root,'logs','mpileup-parts')
    LSF.lsf_run_until_done(to_run_dict,logfile,queue,'-R "select[mem>%s]"' % job_ram, 'mpileup',njobs,MAX_RETRY)

    #logfile = os.path.join(vcf_parts_root,'logs','mpileup-parts-merge')
    #LSF.lsf_run_until_done(merge_subparts_trd,logfile,queue,'-R "select[mem>%s]"' % job_ram, 'mpileup-partsmerge',njobs,MAX_RETRY)
    
    #cmds = unfinished_cmds(to_run_dict)
    #while cmds:
    #    jobids,namedict = LSF.lsf_jobs_submit(cmds,logfile,'normal_serial','-R "select[mem>30000]"', jobname_base='gatk')
    #    LSF.lsf_wait_for_jobs(jobids,logfile,'normal_serial',namedict=namedict)
    #    cmds = unfinished_cmds(to_run_dict)

    cmd = run_safe.safe_script(merge_vcf_parts_cmd(subparts,ref,mpoutvcf,gatk_jar,gatk_ram,tmpdir),mpoutvcf)
    ret = os.system(cmd)
    if ret != 0:
        raise OSError, 'VCF merge failed:\n%s' % cmd



def fqname_from_sample_dict(d,read_subpath='reads_by_individual',readnum=1,read_ext='txt.gz'):
    '''for stampy-fix paths, set (for example):
    read_subpath="reads_by_individual/stampy-fix"

    for read 2, set readnum=2
    '''
    newfq = '%s/%s/%s_lane%s%s/%s_s_%s_%s_sequence%s.%s' % (d['datapath'],read_subpath,d['flowcell'],d['lane'],idxstr_from_idx(d.get('index','')),d['sampleid'],d['lane'],readnum,idxstr_from_idx(d.get('index','')),read_ext)
    return newfq

def bam_from_sample_dict(d,suffix,read_subpath='reads_by_individual'):
    newbam = '%s/%s/%s_lane%s%s/%s_%s' % (d['datapath'],read_subpath,d['flowcell'],d['lane'],idxstr_from_idx(d.get('index','')),d['sampleid'],suffix)
    return newbam

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

        if read.endswith('.merge.fastq.gz'):
            print >> sys.stderr, '%s indicates merged PE reads; skip pairing' % '.merge.fastq.gz'
            unpaired.append(read)
            skip.append(read)
            continue
            
        readpath,readbase = os.path.split(read)
        try:
            if '_s_' in readbase:
                ind,voidstr,lane,suf = re.search('^(.+?)(_s)_(\d)_(.+?)$',readbase).groups()
            else:
                voidstr = ''
                ind,lane,suf = re.search('^(.+?)_(\d)_(.+?)$',readbase).groups()
            fmt = 'legacy'
        except:
            try:
                print >> sys.stderr, 'legacy parse failed for',readbase,'check new'
                ind,suf = re.search('^(.+?)(_index.+?)$',readbase).groups()
                fmt = 'new'
            except:
                try:
                    print >> sys.stderr, 'legacy, new parse failed for',readbase,'check "Sample"'
                    ind,suf = re.search('^(.+?)(_Sample.+?)$',readbase).groups()
                    fmt = 'sample'
                except:
                    print >> sys.stderr, 'parse failed for',readbase
                    raise

        if fmt == 'new':
            readset = set(['1','2'])
            presuf,readnum,postsuf = re.search('(^.+?_L\d{3}_R)(\d)(_.+?)$',suf).groups()
            if readnum in readset:
                r2num = (readset-set(readnum)).pop()
            else:
                print >> sys.stderr, 'read %s not in valid set %s' % (readnum,readset)
                raise ValueError
            matebase = ind+presuf+r2num+postsuf
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
        elif fmt == 'sample':
            readset = set(['1','2'])
            presuf,readnum,postsuf = re.search('(^.+?\.R)(\d)(\..+?)$',suf).groups()
            if readnum in readset:
                r2num = (readset-set(readnum)).pop()
            else:
                print >> sys.stderr, 'read %s not in valid set %s' % (readnum,readset)
                raise ValueError
            matebase = ind+presuf+r2num+postsuf
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
        elif fmt == 'legacy':
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
                matebase = '_'.join([ind+voidstr,lane,'1'+suf[1:]])
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

        else:
            errstr = 'fmt must be "legacy" or "new", is %s' % fmt
            raise ValueError,errstr

    return unpaired,paired

if __name__ == '__main__':

    import argparse

    ds =  ' [%(default)s]'
    #create command line parser
    parser = argparse.ArgumentParser(description='performs stampy mapping')
    parser.add_argument('-v','--vcfname',default=None,help='if vcfname (filename) is supplied, will run GATK with --gatk_argstr and any additional BAM files specified using -a. VCF will be created in OUTROOT'+ds)
    parser.add_argument('-s','--stampy_argstr',default="'--sensitive --substitutionrate=0.02 --maxbasequal=60'",type=eval,help='arguments passed to stampy. Must be single AND double quoted for spaces'+ds)
    parser.add_argument('-g','--gatk_argstr',default="'-out_mode EMIT_ALL_CONFIDENT_SITES -dcov 100'",type=eval,help='arguments passed to GATK (only relevant if --outvcf specified) Must be single AND double quoted for spaces.'+ds)
    parser.add_argument('-gh','--gatkhaplo_argstr',default="'-out_mode EMIT_ALL_CONFIDENT_SITES -dr 50'",type=eval,help='arguments passed to GATKHaplotypeCaller (only relevant if --outvcf specified) Must be single AND double quoted for spaces.'+ds)
    parser.add_argument('-mp','--mpileup_argstr',default="''",type=eval,help='arguments passed to mpileup (only relevant if --outvcf specified) Must be single AND double quoted for spaces, e.g. "\'-r chr3\'"'+ds)
    parser.add_argument('-gr','--gatk_ram',default=8,type=int,help='java VM ram size (in GB).'+ds)
    
    parser.add_argument('-a','--add_bam',default=[],action='append',help='additional BAM files for GATK (only relevant if --outvcf specified) any number of -a arguments may be supplied'+ds)

    parser.add_argument('-r','--reads_per_part',default=100000,type=int,help='target number of reads for each stampy process, used to calculate stampy --processpart.'+ds)
    parser.add_argument('-n','--num_batches',default=100,type=int,help='number of LSF batches to submit.'+ds)
    parser.add_argument('-q','--lsf_queue',default='normal_serial',type=str,help='LSF submission queue.'+ds)
    parser.add_argument('-fbq','--fallback_queue',default='',type=str,help='fallback LSF submission queue; invoked after MAX_RETRY failed LSF reruns on --lsf_queue.'+ds)

    parser.add_argument('-tr','--target_regions',default=None,help='file of sequences (one seq ID per line) to genotype across.'+ds)

    parser.add_argument('--no_merge',action='store_true',help='do not perform bam merge after mapping (regardless of number of individual bams)'+ds)
    parser.add_argument('--fast_merge',action='store_true',help='skip validation in bam merge'+ds)
    
    parser.add_argument('-up','--unpair_reads',action='store_true',help='will NOT attempt to pair read1 and read2 files in input reads'+ds)
    parser.add_argument('--realign',action='store_true',help='perform targeted realignment on individual BAM files'+ds)
    parser.add_argument('--force_realign',action='store_true',help='force realignment (must also specify --realign to use realigned files in genotyping)'+ds)
    parser.add_argument('--reduce_reads',action='store_true',help='invoke individual-level realignment and ReduceReads on each sample BAM'+ds)
    parser.add_argument('--cleanup',action='store_true',help='remove intermediate .sam files created during mapping'+ds)
    parser.add_argument('--debug',action='store_true',help='print commands but do not run anything'+ds)
    parser.add_argument('--bamcheck',action='store_true',help='force checking bam inputs (do not trust previous validation, even if present)'+ds)

    parser.add_argument('--check_donefiles',action='store_true',help='force remove of .done files with missing .bams'+ds)
    
    parser.add_argument('--skip_mpileup',action='store_true',help='do not perform samtools mpileup genotyping, even if -v is supplied (i.e. run only GATK)'+ds)
    parser.add_argument('--skip_haplo',action='store_true',help='do not perform GATK HaplotypeCaller genotyping, even if -v is supplied (i.e. run only GATK)'+ds)

    parser.add_argument('reference_fasta',help='reference for stampy')
    parser.add_argument('outroot',help='directory for logfile and vcf creation')
    parser.add_argument('reads',nargs='+',help='fastq for stampy (qualities should be fastq-33)')
    
    opts = parser.parse_args()

    if opts.reduce_reads:
        opts.realign = False
        opts.force_realign = False
        opts.skip_mpileup = True
        print >> sys.stderr, '--reduce_reads invoked, --skip_mpileup is enabled and the following options are disabled: --realign --force_realign'

    gatkRAM = opts.gatk_ram # gb ram for GATK

    vcfname, stampy_argstr, reference_fasta = opts.vcfname,opts.stampy_argstr,opts.reference_fasta
    njobs = opts.num_batches

    reads = opts.reads

    outroot = opts.outroot
    
    t = reference_fasta
    if not os.path.exists(outroot): os.makedirs(outroot)

    #create indices for reference if not present
    samtools_faidx = reference_fasta+'.fai'
    picard_dict = os.path.splitext(reference_fasta)[0]+'.dict'
    if not os.path.exists(samtools_faidx):
        ret = os.system('samtools faidx %s' % reference_fasta)
        if ret != 0:
            raise OSError, 'unable to create samtools faidx for reference'
    if not os.path.exists(picard_dict):
        ret = os.system('java -jar %s REFERENCE=%s OUTPUT=%s' % (picard_seqdict_jar,reference_fasta,picard_dict))
        if ret != 0:
            raise OSError, 'unable to create picard dict for reference'        


    tb = os.path.splitext(os.path.basename(t))[0]
    bp = stampy_argstr.replace(' ','').replace('=','')

    gidx = t+'.stidx'
    ghash = t+'.sthash'

    if not os.path.exists(gidx):
        #print >> sys.stderr, 'ssh $HOSTNAME "module load %s; stampy.py -G %s %s"' % (stampy_module,t,t)
        os.system('ssh $HOSTNAME "module load %s; stampy.py -G %s %s"' % (stampy_module,t,t))
        #time.sleep(5)
        if not os.path.exists(gidx):
            raise OSError, gidx+' does not exist'
    if not os.path.exists(ghash):
        os.system('ssh $HOSTNAME "module load %s; stampy.py -g %s -H %s"' % (stampy_module,t,t))
        #time.sleep(5)
        if not os.path.exists(ghash):
            raise OSError, ghash+' does not exist'
       
    
    cmds = []
    bams = []
    rg_ref_bams = []
    samparts_by_bam = {}
    cmd_by_sam = {}
    flowcell_by_bamstr = {}

    if opts.unpair_reads:
        print >> sys.stderr, '%s fastq processed (no PE check)' % len(reads)
    else:
        print >> sys.stderr, 'PE check invoked on %s fastq' % len(reads)
        unpaired,paired = find_paired_reads(reads)
        print >> sys.stderr, 'found %s PE; %s SR' % (len(paired),len(unpaired))
        reads = unpaired+paired

    for rfnum,readfile in enumerate(reads):
        print >> sys.stderr, '\n%s / %s : prepare stampy run for %s' % (rfnum+1,len(reads),readfile)
        #throughout, tuple indicates PE
        if isinstance(readfile,tuple):
            readroot,readfilebase = os.path.split(readfile[0])
        else:
            readroot,readfilebase = os.path.split(readfile)

        try:
            rb = re.search('(.+?)_([s\d]|index|Sample)_',readfilebase).groups()[0]
        except:
            print >> sys.stderr, readfilebase
            raise

        fc = os.path.basename(readroot)
        if rb.startswith('pass'): continue

        if isinstance(readfile,tuple):
            samfbase = os.path.join(readroot,'%s_pair-%s_stampy%s' % (rb,tb,bp))
        else:
            samfbase = os.path.join(readroot,'%s-%s_stampy%s' % (rb,tb,bp))

        #sortbam = samf[:-3] + 'sort.bam'
        rg_ref_bam = samfbase + '.rg_refsort.bam'
        rg_ref_bams.append(rg_ref_bam)
        #flowcell_by_bamstr[sortbam] = os.path.dirname(readroot)
        

        #bams.append(sortbam)
        ID = rb+'_'+fc
        SM = rb

        headline = '@RG\tID:%s\tPL:Illumina\tLB:%s\tSM:%s\n' % (ID,SM,SM)
        
        readgroup_arg = 'ID:%s,PL:Illumina,LB:%s,SM:%s' % (ID,SM,SM)

        #print >> sys.stderr, t, readgroup_arg, readfile

        donefile = rg_ref_bam+'.done'
        if os.path.exists(rg_ref_bam) and not os.path.exists(donefile): #not necessarily right -- should be a run-time option
            print >> sys.stderr, 'donefile absent for %s; remove and re-run' % rg_ref_bam
            os.unlink(rg_ref_bam)
            
        if os.path.exists(rg_ref_bam):
            #check bam completeness
            bam_validate = rg_ref_bam+'.valid'
            if os.path.exists(bam_validate) and not opts.bamcheck:
                try:
                    valsize = float(open(bam_validate).read())
                except ValueError:
                    valsize = None
                if valsize and os.path.getsize(rg_ref_bam) == valsize:
                    print >> sys.stderr, 'skip',rg_ref_bam,'(exists; cached validation checked)'
                    continue
            print >> sys.stderr, rg_ref_bam,'found; validate'
            valcmd =  'java -jar /n/home08/brantp/src/picard_svn/trunk/dist/ValidateSamFile.jar INPUT=%s MODE=SUMMARY MAX_OPEN_TEMP_FILES=100 IGNORE=MISMATCH_FLAG_MATE_NEG_STRAND IGNORE=MISMATCH_FLAG_MATE_UNMAPPED IGNORE=MISMATCH_MATE_ALIGNMENT_START IGNORE=MISMATCH_MATE_REF_INDEX 2> /dev/null' % rg_ref_bam
            ret = os.system(valcmd)
            if ret == 0:
                open(bam_validate,'w').write('%s\n' % os.path.getsize(rg_ref_bam))
                print >> sys.stderr, 'skip',rg_ref_bam,'(exists; valid)'
                continue #skip this individual if already mapped under these params
            else:
                print >> sys.stderr, 'failed validation, quit'
                raise OSError

        else:
            if opts.check_donefiles:
                if os.path.exists(donefile):
                    print >> sys.stderr, '.done file for bam %s found, but bam is missing; remove %s ...' % (rg_ref_bam,donefile),
                    ret = os.system('rm -f %s' %  donefile)
                    if ret == 0:
                        print >> sys.stderr, 'DONE'
                    else:
                        raise OSError, 'FAILED'
            
        if isinstance(readfile,tuple):
            r1,r2 = readfile

            readct1 = preprocess_radtag_lane.get_read_count(r1)
            readct2 = preprocess_radtag_lane.get_read_count(r2)
            if readct1 != readct2:
                raise ValueError, 'read counts do not match, abort'

            num_parts = (readct1*2)/opts.reads_per_part
            if num_parts < 1: num_parts = 1
            print >> sys.stderr, 'map in %s part(s)' % num_parts
            
            samparts_by_bam[rg_ref_bam] = []
            for i in xrange(1,num_parts+1):
                sampart = samfbase+'_%05dof%05d.sam' % (i,num_parts)
                samparts_by_bam[rg_ref_bam].append(sampart)
                cmdstr = 'run_safe.py \"module load %s; stampy.py -g %s -h %s  --gatkcigarworkaround --overwrite --readgroup=%s --processpart %s/%s -o %s %s -M %s %s\" %s.done' % (stampy_module,t,t,readgroup_arg,i,num_parts,sampart,stampy_argstr,r1,r2,sampart)
                cmds.append(cmdstr)
                cmd_by_sam[sampart] = cmdstr

        else:
            r1 = readfile

            readct1 = preprocess_radtag_lane.get_read_count(r1)

            num_parts = readct1/opts.reads_per_part
            if num_parts < 1: num_parts = 1
            print >> sys.stderr, 'map in %s part(s)' % num_parts
            
            samparts_by_bam[rg_ref_bam] = []
            for i in xrange(1,num_parts+1):
                sampart = samfbase+'_%05dof%05d.sam' % (i,num_parts)
                samparts_by_bam[rg_ref_bam].append(sampart)
                cmdstr = 'run_safe.py \"module load %s; stampy.py -g %s -h %s  --gatkcigarworkaround --overwrite --readgroup=%s --processpart %s/%s -o %s %s -M %s\" %s.done' % (stampy_module,t,t,readgroup_arg,i,num_parts,sampart,stampy_argstr,r1,sampart)
                cmds.append(cmdstr)
                cmd_by_sam[sampart] = cmdstr

    if opts.debug:
        print 'COMMANDS FOLLOW:\n'+'\n'.join(cmds)
    else:
        logfile = os.path.join(outroot,'lsflog','stampy-%s-%s-log' % (bp,tb))
        LSF.lsf_run_until_done(cmd_by_sam,logfile,'normal_serial','-R "select[mem>20000]"', 'stampy',njobs,MAX_RETRY)
        #if len(cmds) > 0:
        #    
        #    logfile = os.path.join(outroot,'lsflog','stampy-%s-%s-log' % (bp,tb))
        #    
        #    jobids,namedict = LSF.lsf_jobs_submit(cmds,logfile,'normal_serial',jobname_base='stampy',num_batches=njobs)
        #    time.sleep(20)
        #    LSF.lsf_wait_for_jobs(jobids,logfile,namedict=namedict)
        #
        #    cmds = unfinished_cmds(cmd_by_sam)
        #
        #    retries = 0
        #    last_cmds = []
        #    while len(cmds) > 0:
        #        #code to halt execution on recurrent errors
        #        if set(last_cmds) == set(cmds):
        #            if retries > MAX_RETRY:
        #                raise IOError, 'maximum number of retry attempts (%s) exceeded with identical jobs lists.  Check stampy logs for recurrent errors' % MAX_RETRY
        #            else:
        #                retries += 1
        #        last_cmds = cmds
        #        
        #        jobids,namedict = LSF.lsf_jobs_submit(cmds,logfile,'normal_serial',jobname_base='stampy',num_batches=njobs)
        #        time.sleep(20)
        #        LSF.lsf_wait_for_jobs(jobids,logfile,namedict=namedict)
        #
        #        cmds = unfinished_cmds(cmd_by_sam)


    #MERGE SAM PARTS FROM STAMPY
    cmds = []
    mergecmds_by_bam = {}
    for bam,sams in samparts_by_bam.items():
        cmd = 'merge_sams_with_validation.py %s %s' % (bam,' '.join(sams))
        cmds.append(cmd)
        mergecmds_by_bam[bam] = run_safe.safe_script(cmd,bam)

    if opts.debug:
        print 'COMMANDS FOLLOW:\n'+'\n'.join(cmds)
    else:
        logfile = os.path.join(outroot,'lsflog','merge-%s-%s-log' % (bp,tb))
        LSF.lsf_run_until_done(mergecmds_by_bam,logfile,opts.lsf_queue,'-R "select[mem>20000]"','stampy-merge',njobs,MAX_RETRY)
        #if len(cmds) > 0:
        #    import time
        #    logfile = os.path.join(outroot,'lsflog','merge-%s-%s-log' % (bp,tb))
        #    jobids,namedict = LSF.lsf_jobs_submit(cmds,logfile,'normal_serial',jobname_base='stampy-merge',num_batches=njobs)
        #    time.sleep(20)
        #    LSF.lsf_wait_for_jobs(jobids,logfile,namedict=namedict)
        #
        #    cmds = unfinished_cmds(mergecmd_by_bam)
        #    while len(cmds) > 0:
        #        jobids,namedict = LSF.lsf_jobs_submit(cmds,logfile,'normal_serial',jobname_base='stampy-merge')
        #        time.sleep(20)
        #        LSF.lsf_wait_for_jobs(jobids,logfile,namedict=namedict)
        #        cmds = unfinished_cmds(mergecmd_by_bam)
        #
        #    cmds = []

    if opts.cleanup:
        print >> sys.stderr, 'remove %s .sam part files' % len(cmd_by_sam)
        for i,f in enumerate(cmd_by_sam.keys()):
            os.unlink(f)
            os.unlink(f+'.done')
            print >> sys.stderr,'\r%s' % (i+1),

    # MAPPING COMPLETE


    # Genotyping steps follow
    #  IF --reduced_reads, then realign and reduce before merging to single bam
    #  otherwise merge into 1 bam, then realign (if --realign) then perform relevant GATK/samtools/whatever steps
    
    # run reduced reads (and single-sample realignment)
    if opts.reduce_reads:
        to_run_dict = {}
        rr_rg_ref_bams = []
        for bam in rg_ref_bams:
            rr_bam = os.path.splitext(bam)[0] + '.realigned.reduced.bam'
            rr_done = os.path.splitext(bam)[0] + '.realigned.reduced'
            rr_cmd = 'realign_reduce_bam.py %s %s' % (bam,reference_fasta)
            run_safe.add_cmd(to_run_dict,rr_done,rr_cmd,force_write=True)
            rr_rg_ref_bams.append(rr_bam)
        logfile = os.path.join(outroot,'lsflog','realign-reduce-log')
        LSF.lsf_run_until_done(to_run_dict,logfile,opts.lsf_queue,'-R "select[mem>20000]"','realign-reduce',njobs,MAX_RETRY)
        #switch rg_ref_bams to reduced
        orig_rg_ref_bams = rg_ref_bams
        rg_ref_bams = rr_rg_ref_bams

    missing_bams = []
    for b in rg_ref_bams:
        if not os.path.exists(b):
            missing_bams.append(b)
    if len(missing_bams) != 0:
        print >> sys.stderr, 'the following inputs are missing:\n\t%s' % ('\n\t'.join(missing_bams))
        raise ValueError, 'missing bams for merge'

    if vcfname is None:
        print >> sys.stderr, 'alignment complete; subsequent steps only invoked if -v / --vcfname set'
        sys.exit(0)

    # MERGE BAMS IF MORE THAN 100 HERE?
    if vcfname is not None and len(rg_ref_bams) > MERGE_BAMS_ABOVE and not opts.no_merge:
        mergebam = os.path.join(outroot,vcfname+'-all_bam-merged.bam')
        if opts.fast_merge:
            cmd = 'merge_sams_no_validation.py %s %s' % (mergebam,' '.join(rg_ref_bams))
        else:
            cmd = 'merge_sams_with_validation.py %s %s' % (mergebam,' '.join(rg_ref_bams))
        ss = run_safe.safe_script(cmd,mergebam,force_write=True)
        print >> sys.stderr, 'attempt:\n',ss
    
        ret = os.system(ss)
        if ret != 0:
            raise OSError, 'failed to merge bams'
        rg_ref_bams_old = rg_ref_bams
        rg_ref_bams = [mergebam]
        if not os.path.exists(mergebam+'.done'):
            print >> sys.stderr, 'merge invoked, donefile %s not found' % mergebam+'.done'
            raise OSError

    #PERFORM REALIGNMENT IF SELECTED
    if opts.realign:
        #do realignment
        realigned_bams = realign_bams_lsf(rg_ref_bams,reference_fasta,outroot,njobs,min_ind_realign,queue=opts.lsf_queue,gatk_ram=gatkRAM,force_links=opts.force_realign,fallback_queue=opts.fallback_queue)
        
        
    
    #COMPILE BAMS FROM READBASES
    if vcfname is not None:
        if opts.target_regions:
            try:
                include_regions = [l.strip() for l in open(opts.target_regions).readlines()]
                print 'restrict genotyping to %s regions from %s' % (len(include_regions),opts.target_regions)
            except:
                print >> sys.stderr, 'failed to load regions from %s; genotype all sequences' % opts.target_regions
                include_regions = None
        else:
            include_regions = None

        #UnifiedGenotyper
        call_variants_gatk_lsf(rg_ref_bams, reference_fasta, outroot, vcfname, \
                               njobs=opts.num_batches, gatk_program='UnifiedGenotyper', \
                               gatk_args=opts.gatk_argstr, gatk_jar=gatk_jar, gatk_ram=gatkRAM, \
                               tmpdir=None, queue=opts.lsf_queue, job_ram='30000', MAX_RETRY=MAX_RETRY, \
                               include_regions=include_regions,compress_vcf=True,fallback_queue=opts.fallback_queue)
        #HaplotypeCaller
        if not opts.skip_haplo:
            call_variants_gatk_lsf(rg_ref_bams, reference_fasta, outroot, vcfname, \
                                   njobs=opts.num_batches, gatk_program='HaplotypeCaller', \
                                   gatk_args=opts.gatkhaplo_argstr, gatk_jar=gatk2_jar, gatk_ram=gatkRAM, \
                                   tmpdir=None, queue=opts.lsf_queue, job_ram='30000', MAX_RETRY=MAX_RETRY, \
                                   include_regions=include_regions,compress_vcf=True,fallback_queue=opts.fallback_queue)

        if opts.realign:
            call_variants_gatk_lsf(realigned_bams, reference_fasta, outroot, vcfname+'-realign', \
                                   njobs=opts.num_batches, gatk_program='UnifiedGenotyper', \
                                   gatk_args=opts.gatk_argstr, gatk_jar=gatk_jar, gatk_ram=gatkRAM, \
                                   tmpdir=None, queue=opts.lsf_queue, job_ram='30000', MAX_RETRY=MAX_RETRY, \
                                   include_regions=include_regions,compress_vcf=True,fallback_queue=opts.fallback_queue)
            if not opts.skip_haplo:
                call_variants_gatk_lsf(realigned_bams, reference_fasta, outroot, vcfname+'-realign', \
                                       njobs=opts.num_batches, gatk_program='HaplotypeCaller', \
                                       gatk_args=opts.gatkhaplo_argstr, gatk_jar=gatk2_jar, gatk_ram=gatkRAM, \
                                       tmpdir=None, queue=opts.lsf_queue, job_ram='30000', MAX_RETRY=MAX_RETRY, \
                                       include_regions=include_regions,compress_vcf=True,fallback_queue=opts.fallback_queue)
        
        #vcfbasename = vcfname.endswith('.vcf') and vcfname[:-4] or vcfname
        #gatkoutvcf = os.path.join(outroot,'%s-GATK.vcf' % vcfbasename)
        #gatkhaplovcf = os.path.join(outroot,'%s-GATKhaplo.vcf' % vcfbasename)
        #mpileupoutvcf = os.path.join(outroot,'%s-mpileup.vcf' % vcfbasename)
        #istr = ' '.join(['-I '+fn for fn in list(set(rg_ref_bams+opts.add_bam))])
        #bstr = ' '.join([fn for fn in list(set(rg_ref_bams+opts.add_bam))])
        #tmpdir = os.path.join(outroot,'gatk_tmp')
        #gatkcmd = 'java -Xmx%sg -Djava.io.tmpdir=%s -jar  %s -R %s -T UnifiedGenotyper -o %s %s %s' % (gatkRAM,tmpdir,gatk_jar,reference_fasta,gatkoutvcf,opts.gatk_argstr,istr)
        #gatk_sh = os.path.join(outroot,'%s-GATK.sh' % vcfbasename)
        #open(gatk_sh,'w').write(gatkcmd)
        #print >> sys.stderr, 'running GATK:\n%s' % gatkcmd
        #
        #if not opts.debug:
        #    os.system('run_safe.py "bash %s" %s' % (gatk_sh,gatkoutvcf+'.done'))
        #
        #gatkcmd = 'java -Xmx%sg -Djava.io.tmpdir=%s -jar  %s -R %s -T HaplotypeCaller -o %s %s %s' % (gatkRAM,tmpdir,gatk2_jar,reference_fasta,gatkhaplovcf,opts.gatkhaplo_argstr,istr)
        #gatk_sh = os.path.join(outroot,'%s-GATKhaplo.sh' % vcfbasename)
        #open(gatk_sh,'w').write(gatkcmd)
        #print >> sys.stderr, 'running GATK:\n%s' % gatkcmd
        #
        #if not opts.debug:
        #    os.system('run_safe.py "bash %s" %s' % (gatk_sh,gatkhaplovcf+'.done'))        
        
        if not opts.skip_mpileup:
            print >> sys.stderr, 'running mpileup'
            #mpileup
            call_variants_mpileup_lsf(rg_ref_bams, reference_fasta, outroot, vcfname, \
                                      njobs=opts.num_batches, mpileup_args=opts.mpileup_argstr, gatk_jar=gatk_jar, gatk_ram=gatkRAM, \
                                      tmpdir=None, queue=opts.lsf_queue, job_ram='30000', MAX_RETRY=MAX_RETRY, \
                                      include_regions=include_regions)
            if opts.realign:
                call_variants_mpileup_lsf(realigned_bams, reference_fasta, outroot, vcfname+'-realign', \
                                          njobs=opts.num_batches, mpileup_args=opts.mpileup_argstr, gatk_jar=gatk_jar, gatk_ram=gatkRAM, \
                                          tmpdir=None, queue=opts.lsf_queue, job_ram='30000', MAX_RETRY=MAX_RETRY, \
                                          include_regions=include_regions)
