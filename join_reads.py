#!/usr/bin/env python

import Seq, os,sys
from radtag_denovo import preprocess_radtag_lane
from Util import smartopen

def join_pair(r1,r2,num_n=10,qual_n='#'):
    return [r1[0],r1[1]+'N'*num_n+str(Seq.Sequence(r2[1]).rc()),r1[2]+qual_n*num_n+''.join(reversed(r2[2]))]

if __name__ == "__main__":
    f1,f2 = sys.argv[1:]
    fh1 = smartopen(f1)
    fh2 = smartopen(f2)

    rc = preprocess_radtag_lane.get_read_count(f1)
    
    for i in xrange(rc):
        if i % 1000 == 0:
            print >> sys.stderr, '\r%s / %s' % (i,rc),
        r1 = preprocess_radtag_lane.next_read_from_fh(fh1,4)
        r2 = preprocess_radtag_lane.next_read_from_fh(fh2,4)
        print preprocess_radtag_lane.as_fq4_lines(*join_pair(r1,r2))
    print >> sys.stderr, '\ndone'

    
    
