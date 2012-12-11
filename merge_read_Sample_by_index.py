#!/usr/bin/env python

import os,sys,re
from glob import glob
from radtag_denovo import preprocess_radtag_lane

idx_lookup = {'ACAGTG': 5,
 'ACTTGA': 8,
 'ATCACG': 1,
 'CAGATC': 7,
 'CGATGT': 2,
 'CTTGTA': 12,
 'GATCAG': 9,
 'GCCAAT': 6,
 'GGCTAC': 11,
 'TAGCTT': 10,
 'TGACCA': 4,
 'TTAGGC': 3,
 'AGCGAC':	13,
 'AGTTAG':	14,
 'CAGAGC':	15,
 'CGTCTG':	16,
 'CTAATC':	17,
 'CTCTAA':	18,
 'GCTGCT':	19,
 'GTAGTA':	20,
 'GTTAGA':	21,
 'TACACT':	22,
 'TAGTAG':	23,
 'TCAACT':	24,       }

def check_merged_reads(merge_d, source_d):
    results = {}
    for idx_seq,idx in idx_lookup.items():
        for r in [1,2]:
            for l in range(1,9):
                sources = glob('%s/*_%s_L00%s_R%s*.fastq.gz' % (source_d,idx_seq,l,r))
                if len(sources) > 0:
                    source_sum = sum([preprocess_radtag_lane.get_read_count(f) for f in sources])
                    merge_f = '%s/s_%s_%s_sequence_index%s.txt.gz' % (merge_d,l,r,idx)
                    if os.path.exists(merge_f):
                        merge_sum = preprocess_radtag_lane.get_read_count(merge_f)
                        results[merge_f] = source_sum == merge_sum
                    else:
                        results[merge_f] = None
    return results

def outfile_from_fname(fname):
    idx_seq,lane,readnum = re.search('_([ATCG]{6})_L(\d+)_R(\d)_',fname).groups()
    return 's_%s_%s_sequence_index%s.txt' % (int(lane),readnum,idx_lookup[idx_seq])

if __name__ == '__main__':
    d = sys.argv[1]
    fqs = sorted(glob(os.path.join(d,'*.fastq.gz')))
    for fq in fqs:
        cmd = 'gunzip -c %s >> %s' % (fq,outfile_from_fname(fq))
        os.system(cmd)

    os.system('gzip -9 s_*_sequence*.txt')




