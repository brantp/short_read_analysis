#!/usr/bin/env python

from radtag_denovo import preprocess_radtag_lane

import os,sys,numpy

def describe_fastq(filename):
    if preprocess_radtag_lane.smartopen(filename).read(1) == '@':
        lnum = 4
    else:
        lnum = 1
    baseQ = None
    fh = preprocess_radtag_lane.smartopen(filename)
    while baseQ is None:
        n,s,q = preprocess_radtag_lane.next_read_from_fh(fh)
        baseQ = preprocess_radtag_lane.get_baseQ(q)
    fh.close()
    return lnum,baseQ


fname,qualcut = sys.argv[1:]
qualcut = float(qualcut)

lnum,baseQ = describe_fastq(fname)

fh = preprocess_radtag_lane.smartopen(fname)

while 1:
    n,s,qs = preprocess_radtag_lane.next_read_from_fh(fh,lnum)
    if len(n) == 0: break
    qa = numpy.array([ord(c) - baseQ for c in qs])
    if qa.mean() >= qualcut:
        print preprocess_radtag_lane.as_fq_line(n,s,qa,baseQ,lnum),
