#!/usr/bin/env python
#
#split reads into directories based on tags
#takes 2 args: tags file, input file

#def hamm(s1,s2):
#    d = 0
#    for i in range(len(s1)):
#        if s1[i] != s2[i]:
#            d += 1
#    return d

from editdist import distance as hamm
import sys,os,re

low_qual = '[BCDEFHIJ]'
pass_q = 2 #number of low qual bases in tag to allow before passing over read
pass_s = 0 #number of Ns in tag seq to allow before passing over read
pass_hamm = 0 #hamming dist from valid tag to allow before passing over read

tagfile,infile = sys.argv[1:]

tags = [l.strip() for l in open(tagfile).readlines()]

if len(set([len(t) for t in tags])) > 1:
    print >> sys.stderr, 'all tags must be the same length'
    sys.exit()

tl = len(tags[0])

outroot,ext = os.path.splitext(infile)

fhdict = {}
passfh = open(outroot+'_pass'+ext,'w')

for l in open(infile):
    fields = l.split(':')
    s,q = fields[-2:]
    head = fields[:2]
    t = s[:tl]
    tq = q[:tl]
    if len(re.findall(low_qual,tq)) > pass_q or t.count('N') > pass_s or min([hamm(t_this,t) for t_this in tags]) > pass_hamm:
        passfh.write(l)
    else:
        try:
            fhdict[t]['reads'].write(':'.join(head + [s[tl:],q[tl:]]))
            fhdict[t]['tags'].write(':'.join(head + [s[:tl],q[:tl]+'\n']))
        except KeyError:
            fhdict[t] = {'tags':open(outroot+'_'+t+'_tags'+ext,'w'),'reads':open(outroot+'_'+t+'_reads'+ext,'w')}
            fhdict[t]['reads'].write(':'.join(head + [s[tl:],q[tl:]]))
            fhdict[t]['tags'].write(':'.join(head + [s[:tl],q[:tl]+'\n']))

for tfhdict in fhdict.values():
    for v in tfhdict.values():
        v.close()
