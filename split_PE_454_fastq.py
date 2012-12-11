#!/usr/bin/env python

'''given sff_extract output fastq (interleaves forward, reverse, unpaired)

writes _PE_1_sequence.txt.gz , _PE_2_sequence.txt.gz , _sequence.txt.gz prepending <base> (see options)'''

import os,sys

from radtag_denovo import preprocess_radtag_lane 

if __name__ == '__main__':
	infile, outbase, minlen = sys.argv[1:]
	minlen = int(minlen)

	fw_outfile = outbase+'_PE_1_sequence.fastq.gz'
	rv_outfile = outbase+'_PE_2_sequence.fastq.gz'
	sr_outfile = outbase+'_sequence.fastq.gz'

	infh = preprocess_radtag_lane.smartopen(infile)
	fw_out, rv_out, sr_out = [preprocess_radtag_lane.smartopen(f,'w') for f in  [fw_outfile,rv_outfile,sr_outfile]]

	pair = {}
	id = 1
	while id:
		id,s,q = preprocess_radtag_lane.next_read_from_fh(infh,4)
		end = id.split('.')[-1]
		if end not in ('f','r') and len(s) >= minlen:
			sr_out.writelines(preprocess_radtag_lane.as_fq4_lines(id,s,q))
			pair = {}
		else:
			pair[end] = [id,s,q]
			id,s,q = preprocess_radtag_lane.next_read_from_fh(infh,4)
			end = id.split('.')[-1]
			pair[end] = [id,s,q]
			if len(set([v[0].rsplit('.',1)[0] for v in pair.values()])) == 1 and all([len(v[1]) >= minlen for v in pair.values()]):
				try:
					fw_out.writelines(preprocess_radtag_lane.as_fq4_lines(*pair['f']))
					rv_out.writelines(preprocess_radtag_lane.as_fq4_lines(*pair['r']))
				except:
					print >> sys.stderr, pair
			else:
				for this_id, this_s, this_q in pair.values():
					if len(this_s) >= minlen:
						sr_out.writelines(preprocess_radtag_lane.as_fq4_lines(this_id,this_s,this_q))
			pair = {}

