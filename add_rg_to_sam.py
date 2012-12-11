#!/usr/bin/env python

'''takes sam lines on stdin and an id and sm as argv,
adds header line and appends id to each sam line
'''

import sys

ID,SM = sys.argv[1:3]

for l in sys.stdin:
    print l.strip('\n')+'\tRG:Z:'+ID
