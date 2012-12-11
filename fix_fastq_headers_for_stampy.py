#!/usr/bin/env python

import os, sys
from radtag_denovo import preprocess_radtag_lane

infile,outfile = sys.argv[1:]

if not os.path.exists(os.path.dirname(outfile)):
    os.makedirs(os.path.dirname(outfile))

ifh = preprocess_radtag_lane.smartopen(infile)
ofh = preprocess_radtag_lane.smartopen(outfile,'w')

r = preprocess_radtag_lane.next_read_from_fh(ifh)
while r[0]:
    r[0] = '%s %s:%s' % (tuple(r[0].rsplit(':',2)))
    ofh.write(preprocess_radtag_lane.as_fq4_lines(*r))
    r = preprocess_radtag_lane.next_read_from_fh(ifh)

ifh.close()
ofh.close()
