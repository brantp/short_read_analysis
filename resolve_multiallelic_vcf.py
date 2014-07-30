#!/usr/bin/env python

import os,sys
from short_read_analysis import variant_detection

invcf,outvcf = sys.argv[1:]

#this would be where one might tweak the multiallelic resolution parameters,
# see variant_detection docstrings for various multiallelic fuctions
multiallelic_fn = variant_detection.resolve_multiallelic_sd_fn(-0.01,0.5,0.02)

vcf_obj_ma = variant_detection.load_vcf(invcf, multiallelic_sites=multiallelic_fn, write_thresholded_vcf=outvcf, store_only=[])

