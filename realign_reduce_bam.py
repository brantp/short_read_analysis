#!/usr/bin/env python

'''
single bam realignment and reducedreads using GATK
'''

picardRAM = 2
gatk_ram = 4
njobs = 50
#gatk_jar = '/n/home08/brantp/src/GATK-git/dist/GenomeAnalysisTK.jar'
#gatk2_jar = '/n/home08/brantp/src/GenomeAnalysisTK-2.1-8-g5efb575/GenomeAnalysisTK.jar'
#gatk2_jar = '/n/home08/brantp/src/GenomeAnalysisTK-2.2-16-g9f648cb/GenomeAnalysisTK.jar'
gatk2_jar = '/n/home08/brantp/src/GenomeAnalysisTK-2.4-3-g2a7af43/GenomeAnalysisTK.jar'
gatk_jar = gatk2_jar
#picard_jar_root = '/n/home08/brantp/src/picard_svn/trunk/dist'
picard_jar_root = '/n/home08/brantp/src/picard_svn_20130220/trunk/dist'

targetcreator_opts='--maxIntervalSize 5000 -dcov 100'
reducereads_opts = '--minimum_mapping_quality 0'

remove_workdir = False

import os,sys
from map_reads_by_indiv_stampy import partition_reference
import run_safe

bam,ref = sys.argv[1:]

bam_base = os.path.splitext(bam)[0]
workdir = bam_base + '-realign-working/'
if not os.path.exists(workdir):
    os.makedirs(workdir)

intervals_file = bam_base + '.RealignTargets.intervals'
realigned_bam = bam_base + '.realigned.bam'
reduced_bam = bam_base + '.realigned.reduced.bam'

intervals_parts_regions = partition_reference(ref,njobs)
intervals_parts = []

for i,part in enumerate(intervals_parts_regions):
    #reg_str = ' -L '.join(part)
    reg_parts_file = os.path.join(workdir,'part%s.intervals' % (i))
    open(reg_parts_file,'w').writelines([p+'\n' for p in part])
    
    intervals_parts_file = os.path.join(workdir,'part%s.RealignTargets.intervals' % (i))
    rtc_part_cmd = 'java -Xmx%sg -jar %s -T RealignerTargetCreator -I %s -R %s -L %s %s -o %s' % \
                   (gatk_ram,gatk_jar,bam,ref,reg_parts_file,targetcreator_opts,intervals_parts_file)
    rtc_ss = run_safe.safe_script(rtc_part_cmd,intervals_parts_file,force_write=True)
    ret = os.system(rtc_ss)
    if ret == 0:
        print >> sys.stderr, '%s / %s complete' % (i+1,len(intervals_parts_regions))
    else:
        errstr = 'failed on %s' % intervals_parts_file
        raise OSError,errstr
    intervals_parts.append(intervals_parts_file)


cat_cmd = 'cat %s > %s' % (' '.join(intervals_parts),intervals_file)
cat_ss = run_safe.safe_script(cat_cmd,intervals_file,force_write=True)
ret = os.system(cat_ss)
if ret == 0:
    print >> sys.stderr, 'intervals parts concatenation finished'
else:
    errstr = 'interval parts concatenation failed'
    raise OSError,errstr

if remove_workdir and os.path.exists(workdir):
    ret = os.system('rm -rf %s' % workdir)
    if ret != 0:
        raise OSError, 'rm workdir failed'

ir_cmd = 'java -Xmx%sg -jar %s -T IndelRealigner -model USE_SW -I %s -R %s --targetIntervals %s -o %s' % \
         (gatk_ram,gatk_jar,bam,ref,intervals_file,realigned_bam)

ir_ss = run_safe.safe_script(ir_cmd,realigned_bam,force_write=True)
ret = os.system(ir_ss)
if ret == 0:
    print >> sys.stderr, 'IndelRealigner finished'
else:
    errstr = 'IndelRealigner failed'
    raise OSError,errstr


rr_cmd = 'java -Xmx%sg -jar %s -T ReduceReads -I %s -R %s %s -o %s' % \
          (gatk_ram,gatk_jar,realigned_bam,ref,reducereads_opts,reduced_bam)

rr_ss = run_safe.safe_script(rr_cmd,reduced_bam,force_write=True)
print rr_ss
ret = os.system(rr_ss)
if ret == 0:
    print >> sys.stderr, 'ReduceReads finished'
else:
    errstr = 'ReduceReads failed'
    raise OSError,errstr

