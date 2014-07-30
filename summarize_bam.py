#!/usr/bin/env python

'''given a bam file, calculates a dict of basic statistics regarding sequence read counts and mapping


'''

import re,sys,os
from subprocess import Popen,PIPE
from collections import defaultdict

'''
print >> sys.stderr, 'loading category lists'

target_categories = [('agouti',eval(open('/n/hoekstra_lab/sureselect/agouti_keeps.list').read())), \
		     ('mc1r',eval(open('/n/hoekstra_lab/sureselect/mc1r_keeps.list').read())), \
		     ('QD6_genomic',eval(open('/n/hoekstra_lab/sureselect/keeps_all_pol_qd6g200qd12.list').read())), \
		     ('QD20_genomic',eval(open('/n/hoekstra_lab/sureselect/keeps_all_pol_qd20qd20.list').read())), \
		     ('non_target',eval(open('/n/hoekstra_lab/sureselect/ruf_bacs.list').read())), \
		     ('2xCov',eval(open('/n/hoekstra_lab/myselect_design/mean_coverage_2.0.list').read())), \
		     ('3xCov',eval(open('/n/hoekstra_lab/myselect_design/mean_coverage_3.0.list').read())), \
		     ('4xCov',eval(open('/n/hoekstra_lab/myselect_design/mean_coverage_4.0.list').read())), \
		     ('5xCov',eval(open('/n/hoekstra_lab/myselect_design/mean_coverage_5.0.list').read()))]


target_categories = [('agouti',eval(open('/n/hoekstra_lab/sureselect/agouti_keeps.list').read())), \
		     ('mc1r',eval(open('/n/hoekstra_lab/sureselect/mc1r_keeps.list').read())), \
		     ('target',eval(open('/n/hoekstra_lab/sureselect/targets/intersect_keeps.list').read())), \
		     ('non_target',eval(open('/n/hoekstra_lab/sureselect/ruf_bacs.list').read())), \
		     ('2xCov',eval(open('/n/hoekstra_lab/myselect_design/mean_coverage_2.0.list').read())), \
		     ('3xCov',eval(open('/n/hoekstra_lab/myselect_design/mean_coverage_3.0.list').read())), \
		     ('4xCov',eval(open('/n/hoekstra_lab/myselect_design/mean_coverage_4.0.list').read())), \
		     ('5xCov',eval(open('/n/hoekstra_lab/myselect_design/mean_coverage_5.0.list').read()))]

for c,l in target_categories:
	print >> sys.stderr, '%s: %s contigs' % (c,len(l))
'''
		     
genome_size = 3300000000.0 #bases


def parse_samflag(flag):
	flag_head = ['pair','all_aln','unmapped','mate_unmapped','rev','mate_rev','read1','read2','secondary_aln','qc_fail','dup']
	return dict(zip(flag_head,[bool(int(i)) for i in reversed(list(bin(flag)[2:]))]))

if __name__ == '__main__':

	bam = sys.argv[1]
	fc,basename = os.path.split(bam)[-2:]
	indiv = basename.split('_')[0]

	by_contig_fh = open(bam+'.contig.hybselstats','w')
	
	head_dict = dict([(l.strip().split()[1][3:],int(l.strip().split()[2][3:])) for l in Popen('samtools view -H %s | grep "@SQ"' % bam,shell=True,stdout=PIPE).stdout.readlines()])

	print >> sys.stderr, 'count unmapped reads ...',
	unmapped = int(Popen('samtools view -f 4 %s | wc -l' % bam,shell=True,stdout=PIPE).stdout.read().strip())
	print >> sys.stderr, 'count mapped reads ...'
	reads = defaultdict(int)
	map_bases = defaultdict(int)
	aln_bases = defaultdict(int)
	phandle = Popen('samtools view -F 4 %s' % bam, shell=True,stdout=PIPE)
	for l in phandle.stdout:
		fields = l.strip().split()
		reads[fields[2]] += 1
		map_bases[fields[2]] += len(fields[9])
		aln_bases[fields[2]] += sum([int(f[:-1]) for f in re.findall('\d+\w',fields[5]) if f[-1] != 'S'])
		
	mapped = sum(reads.values())
	pct_mapped = (float(mapped)/(mapped+unmapped))*100

	reads_by_category = defaultdict(int)
	map_bases_by_category = defaultdict(int)
	aln_bases_by_category = defaultdict(int)
	target_dict = dict([(l.strip().split()[1][3:],int(l.strip().split()[2][3:])) for l in Popen('samtools view -H %s | grep "@SQ"' % bam,shell=True,stdout=PIPE).stdout.readlines()])
	target_bases_by_category = defaultdict(int)
	
	tot_reads = mapped+unmapped

	'''
	for contig in target_dict:
		by_contig_fh.write('%s\t%s\t%s\t%s\n' % (contig,target_dict[contig],map_bases[contig]/float(target_dict[contig]),aln_bases[contig]/float(target_dict[contig])))
		for category,contig_list in target_categories:
			if contig in contig_list:
				reads_by_category[category] += reads[contig]
				map_bases_by_category[category] += map_bases[contig]
				aln_bases_by_category[category] += aln_bases[contig]
				target_bases_by_category[category] += target_dict[contig]

	'''

	header_row = 'ID\tflowcell\ttotal_reads\tmapped_reads\tpct_mapped'

	'''
	metrics = ['reads','map_bases','aln_bases','target_bases','map_cov','aln_cov','pct_genome','pct_reads','enrichment']
	for category,contig_list in target_categories:
		for metric in metrics:
			header_row += '\t%s_%s' % (category,metric)

	'''
	
	print header_row

	data_row = '%s\t%s\t%s\t%s\t%0.2f' % (indiv,fc,tot_reads,mapped,pct_mapped)

	'''
	for category,contig_list in target_categories:
		fract_cat = target_bases_by_category[category]/genome_size
		fract_reads = float(reads_by_category[category])/tot_reads
		enrichment = fract_reads/fract_cat
		map_coverage = float(map_bases_by_category[category])/target_bases_by_category[category]
		aln_coverage = float(aln_bases_by_category[category])/target_bases_by_category[category]
		data_row += '\t%s\t%s\t%s\t%s\t%0.1f\t%0.1f\t%0.3f\t%0.3f\t%0.1f' % \
			    (reads_by_category[category], \
			     map_bases_by_category[category], \
			     aln_bases_by_category[category], \
			     target_bases_by_category[category], \
			     map_coverage, \
			     aln_coverage, \
			     100*fract_cat, \
			     100*fract_reads, \
			     enrichment)

	'''
	print data_row
