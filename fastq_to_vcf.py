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
overlap_preprocess.py (rtd) is invoked to merge overlapping paired ends
\tand remove adapter sequence readthrough.
preprocess_radtag_lane.py (rtd) is used to demultiplex according to gdoc spreadsheet

for single read data:
preprocess_radtag_lane.py (rtd) is used to demultiplex according to gdoc spreadsheet

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

from rtd import preprocess_radtag_lane,config
from glob import glob

import os,sys,re
import run_safe

try:
    import LSF
except:
    print >> sys.stderr, 'LSF unavailable'
try:
    import SLURM
except:
    print >> sys.stderr, 'SLURM unavailable'
        
def get_ext(fname):
    if fname.endswith('.gz'):
        base,ext = os.path.splitext(fname)
    else:
        base = fname
        ext = ''
    ext = os.path.splitext(base)[1] + ext
    return ext

def sample_fq_from_expected(expected_fq_d):
    fq_by_sample = []
    for k,v in expected_fq_d.items():
        fqs = glob(k)
        if len(fqs) == v:
            fq_by_sample.extend(fqs)
        else:
            errstr = '%s matches to %s (expected %s)' % (len(fqs),k,v)
            raise ValueError,errstr
    return fq_by_sample

multiplex_idx_db = 'DB_multiplex_indices'
tcp_host = 'heroint3'
MAX_RETRY = 3

map_reads_exec = 'map_reads_by_indiv-stampy.py'

if __name__ == '__main__':

    import argparse

    ds =  ' [%(default)s]'
    #create command line parser
    parser = argparse.ArgumentParser(description=argparse_description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-jr','--job_ram',default=6,type=int,help='job ram size (in GB).'+ds)
    
    parser.add_argument('-n','--num_batches',default=100,type=int,help='target number of batches to submit to scheduler.'+ds)
    parser.add_argument('-q','--queue',default='general',type=str,help='submission queue. (Treated as slurm partition if --scheduler=SLURM)'+ds)
    parser.add_argument('-fbq','--fallback_queue',default='',type=str,help='fallback LSF/slurm submission queue; invoked after MAX_RETRY failed LSF reruns on --lsf_queue.'+ds)
    parser.add_argument('-sched','--scheduler',default='slurm',type=str,help='Scheduler to submit jobs to.  Current support for "lsf" and "slurm"'+ds)
    parser.add_argument('-mjd','--max_job_duration',default=1440,type=int,help='slurm job max duration (in minutes)'+ds)

    parser.add_argument('--force_db_id',action='store_true',help='force mouse database ids for individuals \n(replacing legacy sampleid from DB_library_data)'+ds)

    parser.add_argument('-v','--vcfname',default=None,help='Overrride default vcfname for map_reads_by_indiv-stampy.py \n("-" delimited list of PROJECTS if None)'+ds)

    #argument strings for map_read_by_indiv-stampy.py
    parser.add_argument('-s','--stampy_argstr',default="'--sensitive --substitutionrate=0.02 --maxbasequal=60'",type=eval, \
                        help='arguments passed to stampy. \nMust be single AND double quoted for spaces'+ds)
    parser.add_argument('-g','--gatk_argstr',default="'-out_mode EMIT_ALL_CONFIDENT_SITES -dcov 200 -glm BOTH'",type=eval, \
                        help='arguments passed to GATK UnifiedGenotyper. \nMust be single AND double quoted for spaces.'+ds)
    parser.add_argument('-gh','--gatkhaplo_argstr',default="'-out_mode EMIT_ALL_CONFIDENT_SITES -dr 50'",type=eval, \
                        help='arguments passed to GATK HaplotypeCaller. \nMust be single AND double quoted for spaces.'+ds)
    parser.add_argument('-mp','--mpileup_argstr',default="''",type=eval, \
                        help='arguments passed to mpileup. \nMust be single AND double quoted for spaces, e.g. "\'-r chr3\'"'+ds)

    parser.add_argument('-mr','--mapreads_argstr',default="''",type=eval, \
                        help='additional arguments for map_reads_by_indiv-stampy.py. \nMust be single AND double quoted for spaces, e.g. "\'--cleanup --fast_merge --reduce_reads\'"'+ds)
    
    parser.add_argument('reference_fasta',help='reference for stampy')
    parser.add_argument('outroot',help='directory for logfile and vcf creation')
    parser.add_argument('projects',nargs='+',help='project names from DB_library_data to include in run')

    opts = parser.parse_args()

    if opts.vcfname is None:
        vcfname = '-'.join(opts.projects)
    else:
        vcfname = opts.vcfname

    index_lookup = preprocess_radtag_lane.no_net_get_table_as_dict(multiplex_idx_db,tcp_host)
    td = preprocess_radtag_lane.no_net_get_table_as_dict(config.LIBRARY_DATA,tcp_host)

    td = [d for d in td if d.get('project',None) in opts.projects and d.has_key('datapath')]

    print >> sys.stderr, '%s individual records found for projects %s' % (len(td),opts.projects)

    preprocess_targets = []
    expected_fq_d = {}

    if opts.force_db_id:
        transtable,failures = preprocess_radtag_lane.get_legacy_to_DB_lookup(td)
    
    for d in td: #UPDATE FOR DB ID LOOKUP
        if opts.force_db_id:
            if d['sampleid'] in transtable:
                ind = transtable[d['sampleid']]
            else:
                print >> sys.stderr, '%s not in transtable' % d['sampleid'] 
        else:
            ind = d['sampleid']
            
        r1,r2 = preprocess_radtag_lane.fq_path_from_db_dict(d,index_lookup)
        if r1:
            print >> sys.stderr, r1,r2,d,index_lookup
            preprocess_targets.append(((r1,r2),(d['flowcell'],d['lane'],d.get('index',None),d.get('cutsite','AATTC'))))
            if r2:
                glob_key = os.path.join(d['datapath'], \
                                        'reads_by_individual', \
                                        '%s_lane%s_index%s_trim' % (d['flowcell'],d['lane'],d.get('index',None)), \
                                        '%s_*%s' % (ind,get_ext(r1)) )
                expected_fq_d[glob_key] = 2
                glob_key = os.path.join(d['datapath'], \
                                        'reads_by_individual', \
                                        '%s_lane%s_index%s_merge' % (d['flowcell'],d['lane'],d.get('index',None)), \
                                        '%s_*%s' % (ind,get_ext(r1)) )
                expected_fq_d[glob_key] = 1
            else:
                raise ValueError, 'single reads unsupported'
                #expected_fq_d[glob_key] = 1
        else:
            errstr = 'no fastq for %s' % d
            raise ValueError, errstr

    preprocess_targets = list(set(preprocess_targets))

    ol_to_run_dict = {}
    pp_to_run_dict = {}

    for (r1,r2),(fc,l,idx,cs) in preprocess_targets:
        if r2 is None: #single read; preprocess only
            cmd = 'preprocess_radtag_lane.py -w -u -s %s -fc %s -l %s -idx %s %s %s' % (cs,fc,l,idx,opts.force_db_id and '--force_db_id' or '',r1)
            ss_base = os.path.join(os.path.dirname(r1),'sr_preprocess_lane%s_index%s_DBID%s' % (l,idx,opts.force_db_id))
            run_safe.add_cmd(pp_to_run_dict, ss_base, cmd, force_write=True)
        else:
            cmd = 'overlap_preprocess.py -fc %s -l %s -idx %s -pp "-w -u -s %s %s" %s %s' % (fc,l,idx,cs,opts.force_db_id and '--force_db_id' or '',r1,r2)
            ss_base = os.path.join(os.path.dirname(r1),'ol_preprocess_lane%s_index%s_DBID%s' % (l,idx,opts.force_db_id))
            run_safe.add_cmd(ol_to_run_dict, ss_base, cmd, force_write=True)

    print pp_to_run_dict
    print
    print ol_to_run_dict

    to_run_dict = {}
    to_run_dict.update(pp_to_run_dict)
    to_run_dict.update(ol_to_run_dict)

    jobname_base = 'preprocess'
    logbase = os.path.join(opts.outroot,'slurmlog','preprocess')
    print >> sys.stderr, 'run %s logs in %s' % (jobname_base,logbase)
    SLURM.run_until_done(to_run_dict,jobname_base,logbase,opts.max_job_duration,(opts.job_ram+1)*1024,opts.num_batches,opts.queue,MAX_RETRY=MAX_RETRY)

    #collect individual fastq/fastq pairs
    #(LATER: GET INDIVIDUAL FILES BY DB LOOKUP; REQUIRES HANDLING MOUSE DB ID LOOKUP IF SET)
    fq_to_run = sample_fq_from_expected(expected_fq_d)
    #print fq_to_run

    map_reads_cmd = map_reads_exec + ' -gr %s -n %s -q %s -sched %s -mjd %s -v %s -s \'"%s"\' -g \'"%s"\' -gh \'"%s"\' -mp \'"%s"\' %s %s %s ' % \
                    (opts.job_ram, \
                     opts.num_batches, \
                     opts.queue, \
                     opts.scheduler, \
                     opts.max_job_duration, \
                     vcfname, \
                     opts.stampy_argstr, \
                     opts.gatk_argstr, \
                     opts.gatkhaplo_argstr, \
                     opts.mpileup_argstr, \
                     opts.mapreads_argstr, \
                     opts.reference_fasta, \
                     opts.outroot)

    map_reads_cmd += ' '.join(fq_to_run)
    
    print >> sys.stderr, 'run map_reads in %s' % (os.path.join(opts.outroot,vcfname))
    map_reads_ss = run_safe.safe_script(map_reads_cmd,os.path.join(opts.outroot,vcfname),force_write=True)
    ret = os.system(map_reads_ss)
