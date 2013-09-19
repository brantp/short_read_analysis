#!/usr/bin/env python
'''Control script for execution starting with index/lane level fastq files (e.g. received from RTA/OLB)

See argparse description or try fastq_to_vcf.py -h for more details
'''

argparse_description = '''--------------------------------------------------------------------------------

Pipeline for handling raw GA/hiseq read data to generate genotype vcf(s)

--------------------------------------------------------------------------------

-INPUT
fastq filenames must be in either "legacy" or "standard" format as follows.

legacy:
<datapath>/s_<lane>_<read>_sequence[_index<index>].txt.gz
standard:
<datapath>/Sample_lane<lane>[<optional>]_<index>.R<read>.fastq.gz

where:
<datapath> <lane> <index> are as per fields of same name in google drive
\t(see config.py in your rtd clone for google spreadsheet name)
<read> is one of {1,2} indicating forward or reverse sequencing read
<optional> and anything appearing in [] are not required
\t(e.g. an index is not required in legacy filenames)

--------------------------------------------------------------------------------

-STEPS

--PREPROCESS (overlap/trim/demultiplex)
for paired end data:
overlap_preprocess.py is invoked to merge overlapping paired ends
\tand remove adapter sequence readthrough.
preprocess_radtag_lane.py is used to demultiplex according to gdoc spreadsheet

for single read data:
preprocess_radtag_lane.py is used to demultiplex according to gdoc spreadsheet

--ALIGNMENT
All demultplexed reads are submitted to map_reads_by_indiv-stampy.py
\tstampy (currently no BWA premap) maps in defined read-count subparts
\tpicard merge_sam and validate are invoked to produce one BAM per sample
\t(defined as flowcell/lane/index/barcode[/merge,trim])

--GATK PROCESSING
Optional GATK steps pre-genotyping (recalibration, realignment, reducereads)
\thandled as specified in map_read_by_indiv-stampy.py

--GENOTYPE CALLING
A single multisample genotyping run (parallelized as specified) is performed
\tfor each genotyper requested in map_read_by_indiv-stampy.py invocation

--------------------------------------------------------------------------------

-OUTPUT
lots.  But a single vcf per genotyper is generated in specified outroot
(more documentation of outputs to come)

--------------------------------------------------------------------------------
'''

if __name__ == '__main__':

    import argparse

    ds =  ' [%(default)s]'
    #create command line parser
    parser = argparse.ArgumentParser(description=argparse_description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-gr','--gatk_ram',default=8,type=int,help='java VM ram size (in GB).'+ds)
    
    parser.add_argument('-n','--num_batches',default=100,type=int,help='number of LSF batches to submit.'+ds)
    parser.add_argument('-q','--lsf_queue',default='normal_serial',type=str,help='LSF submission queue. (Treated as slurm partition if --scheduler=SLURM)'+ds)
    parser.add_argument('-fbq','--fallback_queue',default='',type=str,help='fallback LSF/slurm submission queue; invoked after MAX_RETRY failed LSF reruns on --lsf_queue.'+ds)
    parser.add_argument('-sched','--scheduler',default='lsf',type=str,help='Scheduler to submit jobs to.  Current support for "lsf" and "slurm"'+ds)
    parser.add_argument('-mjd','--max_job_duration',default=2880,type=int,help='slurm job max duration (in minutes)'+ds)

    parser.add_argument('reference_fasta',help='reference for stampy')
    parser.add_argument('outroot',help='directory for logfile and vcf creation')
    parser.add_argument('reads',nargs='+',help='fastq for stampy (qualities should be fastq-33)')
    
    opts = parser.parse_args()
