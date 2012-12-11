#!/usr/bin/env python

import sys, Util
from short_read_analysis import extract_genotypes_from_mclgr
from collections import defaultdict


def increment_lg(maploci,increment):
    for loc,(lg,mp) in maploci.items():
        maploci[loc] = (lg+increment,mp)

    return maploci

def merge_maps(maps):
    all_genotypes = defaultdict(dict)
    all_maploci = {}
    increment = 0
    
    for m in maps:
        if ',' in m:
            mapf,mIDf = m.split(',')
        else:
            mapf = m
            mIDf = None
        maploci,genotypes = extract_genotypes_from_mclgr.load_cross_radtag_genotypes(mapf,mIDf)
        #print >> sys.stderr, m,'\n',[(k,len(v)) for k,v in genotypes.items()]


        all_maploci.update(increment_lg(maploci,increment))
        for k,v in genotypes.items():
            all_genotypes[k].update(v)
        increment = max([v[0] for v in all_maploci.values()])

    return all_maploci,all_genotypes

if __name__ == '__main__':

    out_to = sys.argv[1]
    if out_to == '-':
        outfh = sys.stdout
    else:
        outfh = open(out_to,'w')
    maps = sys.argv[2:]
    all_maploci,all_genotypes = merge_maps(maps)

    extract_genotypes_from_mclgr.output_cross_radtag_genotypes(all_maploci,all_genotypes,outfh)
    
