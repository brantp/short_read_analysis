#!/usr/bin/env python

import sys
strip_tag = True

for l in open(sys.argv[1]):
    head,s,q = l.rsplit(':',2)
    head = head.rsplit(':',2)[0]
    print '@%s\n%s\n+\n%s' % (head,s,q),
