#!/usr/bin/env python

from radtag_denovo import preprocess_radtag_lane
import sys

print >> sys.stderr, 'get data'
tbl = preprocess_radtag_lane.get_table_as_dict('DB_library_data',suppress_fc_check=True)
print >> sys.stderr, 'done'

runs = sorted(list(set([(d.get('flowcell',''),d.get('lane',''),d.get('index',''),d.get('adaptersversion','')) for d in tbl if d.get('lib','') == 'RAD'])))

print runs
