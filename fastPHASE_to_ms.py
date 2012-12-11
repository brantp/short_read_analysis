#!/usr/bin/env python
'''
fastPHASE_to_ms.py fastphase.output.ext fastphase_dir posfile_dir ms_outdir

given extension on fastphase files:

dir/phased/conditions_agouti_bac.inp.phased
dir/pos/conditions_agouti_bac.pos    <--- created by variant_detection.write_fastPHASE_genotypes

run with:
fastPHASE_to_ms.py inp.phased dir/phased dir/pos dir/ms 

will deposit .ms files in dir/ms

'''

import os,sys
from glob import glob


ext, fastphase_dir, posfile_dir, ms_outdir, replace_imputed = sys.argv[1:]

if not os.path.exists(ms_outdir): os.makedirs(ms_outdir)

drop_imputed = replace_imputed != 'False'

globstr = os.path.join(fastphase_dir,'*%s' % ext)
phased = sorted(glob(os.path.join(fastphase_dir,'*%s' % ext)))

for ff in phased:
    pf = os.path.join(posfile_dir,os.path.basename(ff)[:-1*len(ext)]+'.pos')
    mf = os.path.join(ms_outdir,os.path.basename(ff)[:-1*len(ext)]+'.ms')
    print >> sys.stderr, ff,pf,mf
    posline = open(pf).readline()
    in_gt = False
    ofh = open(mf,'w')
    header_done = False
    for l in open(ff):
        if l.startswith('END GENOTYPES'):
            print >> sys.stderr, 'done'
            in_gt = False
        if in_gt:
            if not l.strip().startswith('#'):
                gt = l.strip().split()

                if not header_done:
                    if len(posline.strip().split()) != len(gt):
                        raise ValueError, 'lengths of .pos (%s) and phased gt (%s) inconsistent, abort!' % (len(posline.strip().split()),len(gt))
                    ofh.write('//\n')
                    ofh.write('segsites: %s\npositions: %s' % (len(gt),posline))
                    header_done = True
                    
                if drop_imputed:
                    for i,g in enumerate(gt):
                        if g.startswith('['):
                            gt[i] = replace_imputed
                else:
                    gt = [g.strip('[]') for g in gt]

                ofh.write(''.join(gt))
                ofh.write('\n')

        if l.startswith('BEGIN GENOTYPES'):
            print >> sys.stderr, 'start'
            in_gt = True
    ofh.close()
        
