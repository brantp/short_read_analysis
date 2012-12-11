#!/usr/bin/env python
# 
# takes a tags file for preexisting tags, length of desired new tags, and a number of tags to generate

import sys
from editdist import distance
from numpy import array

def skewscore(s,li):
    sksc = 0
    for i,c in enumerate(s):
        counts = {}.fromkeys(['A','T','C','G'],0)
        for seq in li:
            counts[seq[i]] += 1
        if s[i] == sorted(counts.items(),key=lambda x:x[1],reverse=True)[0][0]:
            sksc += 1
    return sksc

oldtags,tl, mindist, want = sys.argv[1:]

if oldtags == 'None':
    tags = []
else:
    tags = open(oldtags).readlines()

mindist = int(mindist)
maxrep = 2
tl = int(tl)
bases = list("ACGT")

bstr = '+'.join(['bases[%s]']*tl) % tuple([chr(itr) for itr in range(97,97+tl)])
forstr = ' '.join(['for %s in range(4)']*tl) % tuple([chr(itr) for itr in range(97,97+tl)])

candidates = eval('[%s %s if not any([(%s).count(base) > %s for base in bases])]' % (bstr,forstr,bstr,maxrep))


newtags = []
for c in candidates:
    if len(tags+newtags) == 0:
        newtags.append(c)
    elif min([distance(s,c) for s in tags+newtags]) >= mindist:
        newtags.append(c)

while len(newtags) > int(want):
    skews = array([skewscore(s,newtags+tags) for s in newtags])
    newtags.pop(skews.argmax())

print '\n'.join( newtags)
