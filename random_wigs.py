#!/usr/bin/env python

outdir_suf = '-permute'

import sys,os,random

def permute_pheno(pfile,outfile):
    plist = open(pfile).read().strip().split()
    random.shuffle(plist)
    open(outfile,'w').write(' '.join(plist) + '\n')

if __name__ == "__main__":
    geno,pheno,suf = sys.argv[1:]
    basedir = os.path.splitext(pheno)[0]+outdir_suf
    if not os.path.exists(basedir): os.makedirs(basedir)
    outbase = os.path.join(basedir,'output-%s-' % suf)
    pp = os.path.join(basedir,'pheno-%s.txt' % suf)
    permute_pheno(pheno,pp)
    ret = os.system('wigs simevo -g %s -s %s -f %s -x 40000 -b 500 -t 100 -d 1' % (geno,pp,outbase))
    if ret != 0:
        raise OSError
    else:
        print >> sys.stderr, 'DONE'
