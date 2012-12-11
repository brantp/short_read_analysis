#!/usr/bin/env python
import os,sys,Util,LSF,numpy
from short_read_analysis import extract_genotypes_from_mclgr
from collections import defaultdict

def run_parallel_blasts(clids,mapname,gr,tab,uniqued,blastdb,grname=None,nbatches = 10):

    if grname is None:
        grname = os.path.basename(gr)

    batchlen = len(clids) / nbatches
    clid_batches = [clids[i:i+batchlen] for i in range(0,len(clids),batchlen)]

    cmds = []
    blastouts = []
    for cb in clid_batches:
        seqfile = os.path.join(outroot,'%s_%s-%s_%sclust_%s.fa' % (os.path.splitext(mapname)[0],cb[0],cb[-1],len(cb),os.path.splitext(grname)[0]))
        blast_out = os.path.join(outroot,'%s-%s.blast' % tuple([os.path.splitext(os.path.basename(f))[0] for f in [seqfile,blastdb]]))
        blastouts.append(blast_out)
        if os.path.exists(blast_out): continue
        clidlist = '\n'.join(cb)
        cmds.append('echo -e "%s" | blast_from_rad_clusters.py %s %s %s %s %s %s' % (clidlist,gr,tab,uniqued,seqfile,blastdb,blast_out))

    logfile = os.path.join(outroot,'rad-clust-blast-log')
    jobids,namedict = LSF.lsf_jobs_submit(cmds,logfile,'normal_serial',bsub_flags='-R "select[mem>30000]"', jobname_base='radclustblast',num_batches=nbatches)
    LSF.lsf_wait_for_jobs(jobids,logfile,namedict=namedict,restart_z=24)

    return blastouts

def get_chrom_size(blastdb):
    try:
        chrom_size = sorted([(k,int(v)) for k,v in [l.strip().split()[:2] for l in open(blastdb+".fai")]],key = lambda x: x[1], reverse=True)
    except:
        os.system('samtools faidx '+blastdb)
        chrom_size = sorted([(k,int(v)) for k,v in [l.strip().split()[:2] for l in open(blastdb+".fai")]],key = lambda x: x[1], reverse=True)

    return chrom_size

def calc_chrom_step(chrom_size,maploci):

    total_size = sum([v for k,v in chrom_size])
    return total_size/len(maploci)
   
def get_hits_in_ref_window(hits_by_site,chrom,start,end):

    hits_in_win = []
    for k in hits_by_site.keys():
        site,numread = k.split('-')
        num,read = numread[:-1],numread[-1]
        for hitchr,(hitpos,hitp) in hits_by_site[k].items():
            if hitchr == chrom and start <= hitpos < end:
                hits_in_win.append((site,num,read,-1*numpy.log10(hitp)))
    return hits_in_win

def populate_synteny_matrix(maploci,hits_by_site,chrom_size,chrom_step):

    lol = []
    rows = []
    cols = []

    # set up rows
    lastchr = None
    lastpos = None
    for site,(lg,gp) in sorted(maploci.items(), key = lambda x:x[1]):
        if lg != lastchr:
            rows.append(('LG_%s' % lg, None))
        rows.append(('%s:%0.2f_%s' % (lg,gp,site),site))
        lastchr,lastpos = lg,gp

    # set up columns
    for chrom,size in chrom_size:
        cols.append(('Chr_%s' % chrom,None))
        for i in range(0,size,chrom_step):
            cols.append(('%s:%d' % (chrom,i),(chrom,i)))    

    for rlab,site in rows:
        if site is None:
            lol.append(numpy.array([-1]*len(cols)))
        else:
            rval = []
            for clab,cpos in cols:
                if cpos is None:
                    rval.append(-1)
                    lastpos = 0
                else:
                    chrom,pos = cpos
                    # gather hits
                    this_hbs = dict([(k,v) for k,v in hits_by_site.items() if k.startswith(site)])
                    winhits = get_hits_in_ref_window(this_hbs,chrom,lastpos,pos)
                    # tally hits
                    hitsums = defaultdict(float)
                    for s,num,read,logp in winhits:
                        hitsums[num] += logp
                    winscore = max(hitsums.values()+[0])
                    rval.append(winscore)
                    lastpos = pos
                    #if winscore: print >> sys.stderr, 'for %s take %s' % (winhits,winscore)
            lol.append(numpy.array(rval))
                    

    return numpy.array(lol),[k for k,v in rows],[k for k,v in cols]
    

def plot_synteny_matrix(mat,rows,cols,fig=1,fontsize='xx-small'):
    import pylab
    pylab.figure(fig)
    pylab.clf()
    pylab.pcolormesh(mat)

    ax = pylab.figure(fig).axes[0]
    ax.set_xticks(range(len(cols)))
    ax.set_yticks(range(len(rows)))
    ax.set_xticklabels(cols)
    ax.set_yticklabels(rows)
    pylab.matplotlib.pyplot.xticks(rotation=90)
    pylab.matplotlib.pyplot.yticks(fontsize=fontsize)
    pylab.matplotlib.pyplot.xticks(fontsize=fontsize)
    

if __name__ == "__main__":

    mapf,id_header,gr,tab,uniqued,blastdb = sys.argv[1:7]
    outroot,mapname = os.path.split(mapf)
    
    
    maploci,genotypes = extract_genotypes_from_mclgr.load_cross_radtag_genotypes(mapf,'skip',id_header=id_header)

    clids = [k for k,v in sorted(maploci.items(),key = lambda x: x[1])]

    blastouts = run_parallel_blasts(clids,mapname,gr,tab,uniqued,blastdb)


    hits_by_site = defaultdict(dict)
    for f in blastouts:
        for l in open(f):
            if l.startswith("#"): continue
            fields = l.strip().split()
            hits_by_site[fields[0]][fields[1]] = (int(fields[8]),float(fields[10]))
        
    chrom_size = get_chrom_size(blastdb)
    chrom_step = calc_chrom_step(chrom_size,maploci)

    mat,rows,cols = populate_synteny_matrix(maploci,hits_by_site,chrom_size,chrom_step)
    
    '''
    for site in [k for k,v in maploci.items() if v[0] == 3]:
    print site, [max([(mc,hit[0],-1*log10(hit[1])) for mc,hit in v.items()],key = lambda x: x[2]) for k,v in chromhits.items() if k.startswith(site)]

    
    '''
