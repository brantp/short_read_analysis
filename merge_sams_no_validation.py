#!/usr/bin/env python

import os,sys,run_safe

bam,sams = sys.argv[1],sys.argv[2:]

#picard_root = '/n/home08/brantp/src/picard_svn/trunk/dist/'
picard_root = '/n/home08/brantp/src/picard_svn_20130220/trunk/dist/'
picardRAM = 2
max_temp = 1000
max_records = 1500000
MAX_PER_RUN = 100
RM_SAMS = True #always overridden to True for sam/bams created as intermediates in large merge sets

#revolting hack
if len(sams) > MAX_PER_RUN:
    sams1 = sams[:len(sams)/2]
    sams2 = sams[len(sams)/2:]
    bam1 = bam+'-1.bam'
    bam2 = bam+'-2.bam'
    for b,s in [(bam1,sams1),(bam2,sams2)]:
        print >> sys.stderr, '\npre-merge %s (%s parts)\n' % (b,len(s))
        cmd = 'merge_sams_no_validation.py %s %s' % (b,' '.join(s))
        ss = run_safe.safe_script(cmd,b,force_write=True)
        print >> sys.stderr, ss
        ret = os.system(ss)
        if ret != 0:
            raise OSError, 'pre-merge failed'
    sams = [bam1,bam2]
    RM_SAMS = True
    
mergecmd = 'java -Xmx%sg -jar %sMergeSamFiles.jar INPUT=%s OUTPUT=%s MERGE_SEQUENCE_DICTIONARIES=true VALIDATION_STRINGENCY=LENIENT; samtools index %s' % (picardRAM,picard_root,' INPUT='.join(sams), bam, bam)
ret = os.system(mergecmd)
if ret == 0:
    print >> sys.stderr, '\nmerge complete\n'
else:
    print >> sys.stderr, '\nfailed:\n',mergecmd
    raise OSError, 'merge failed for %s' % bam

#valcmd =  'java -Xmx%sg -jar %sValidateSamFile.jar INPUT=%s MODE=SUMMARY MAX_OPEN_TEMP_FILES=%s MAX_RECORDS_IN_RAM=%s IGNORE=MISMATCH_FLAG_MATE_NEG_STRAND IGNORE=MISMATCH_FLAG_MATE_UNMAPPED IGNORE=MISMATCH_MATE_ALIGNMENT_START IGNORE=MISMATCH_MATE_REF_INDEX IGNORE=INVALID_INSERT_SIZE' % (picardRAM,picard_root,bam,max_temp,max_records)

#print >> sys.stderr, valcmd
#ret = os.system(valcmd)
#if ret == 0:
#    open(bam+'.valid','w').write('%s\n' % os.path.getsize(bam))

if ret == 0:
    if RM_SAMS:
        for sam in sams:
            try:
                rm_ret = os.unlink(sam)
                rm_ret = os.unlink(sam+'.bai')
                rm_ret = os.unlink(sam+'.done')
                rm_ret = os.unlink(sam+'.sh')
                rm_ret = os.unlink(sam+'.valid')
            except:
                pass
else:
    raise OSError, 'validation failed for %s' % bam
