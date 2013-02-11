#!/usr/bin/env python

import os,sys

bam,sams = sys.argv[1],sys.argv[2:]

picard_root = '/n/home08/brantp/src/picard_svn/trunk/dist/'
picardRAM = 4
max_temp = 500
max_records = 5000000

mergecmd = 'java -Xmx%sg -jar %sMergeSamFiles.jar INPUT=%s OUTPUT=%s MERGE_SEQUENCE_DICTIONARIES=true VALIDATION_STRINGENCY=LENIENT; samtools index %s' % (picardRAM,picard_root,' INPUT='.join(sams), bam, bam)
ret = os.system(mergecmd)
if ret == 0:
    print >> sys.stderr, '\nmerge complete\n'
else:
    print >> sys.stderr, '\nfailed:\n',mergecmd
    raise OSError, 'merge failed for %s' % bam

valcmd =  'java -Xmx%sg -jar %sValidateSamFile.jar INPUT=%s MODE=SUMMARY MAX_OPEN_TEMP_FILES=%s MAX_RECORDS_IN_RAM=%s IGNORE=MISMATCH_FLAG_MATE_NEG_STRAND IGNORE=MISMATCH_FLAG_MATE_UNMAPPED IGNORE=MISMATCH_MATE_ALIGNMENT_START IGNORE=MISMATCH_MATE_REF_INDEX IGNORE=INVALID_INSERT_SIZE 2> /dev/null' % (picardRAM,picard_root,bam,max_temp,max_records)

print >> sys.stderr, valcmd
ret = os.system(valcmd)
if ret == 0:
    open(bam+'.valid','w').write('%s\n' % os.path.getsize(bam))
else:
    raise OSError, 'validation failed for %s' % bam
