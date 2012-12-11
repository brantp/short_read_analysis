#!/usr/bin/env python
'''
simulate_rad_recovery.py genome.fa fract_total num_reads num_indiv site1 site2 sizing_mean sizing_std outroot

simulates a cut from a hypothetical genome of which genome.fa is fract_total
( i.e. 1 megabase of BAC sequence from a 3 gigabase genome has fract_total = 0.00033 )

generates num_indiv draws from a normal with mean:
num_reads / num_indiv
and stdev = mean/3

using recognition sites site1 and site2 (use quotes for regex representing ambiguity, i.e. ApoI (RAATTY) is "[AG]AATT[CT]")
Not yet implemented: for only one enzyme, supply site2 as "None"

and simulating size selection with mean sizing_mean, stdev sizing_std

writes:



'''

import Seq, Util
import os, sys, re, random
import pylab, numpy
import cPickle

from collections import defaultdict

def cut_fa_obj(fa,site1,site2=None):
    frags = []
    for k,s in fa.items():
        print >> sys.stderr, k,len(s)
        try:
            if site2 is not None:
                frags.extend([len(rs) for rs in re.sub(site2,'2|2',re.sub(site1,'1|1',str(s).upper())).split('|')  if set([rs[0],rs[-1]]) == set(['1', '2'])])
            else:
                frags.extend([len(rs) for rs in re.sub(site1,'1|1',str(s).upper()).split('|')  if set([rs[0],rs[-1]]) == set(['1'])])
        except:
            print >> sys.stderr, 'none'

    return frags
    

def cut_fa_obj_re(fa,site1,site2=None):
    frags = []
    frags2 = []
    for k,s in fa.items():
        print >> sys.stderr, k,len(s)
        try:
            if site2 is not None:
                frags.extend([len(rs) for rs in reduce(lambda x,y:x+y,[re.split(site2,cs) for cs in re.split(site1,str(s)) if site2 in cs])])
                frags2.extend([len(rs) for rs in reduce(lambda x,y:x+y,[cs.split(site2) for cs in str(s).split(site1) if site2 in cs])])
            else:
                frags.extend(re.split(site1,str(s)))
        except:
            print >> sys.stderr, 'none'

    return frags,frags2

def simulate_coverage(frags,mean,stdev,num,lower_bound,upper_bound):
    ht,bt = numpy.histogram(frags,bins = numpy.arange(lower_bound,upper_bound))
    h,b = numpy.histogram(numpy.random.normal(mean,stdev,size=num),bins=bt)
    coverage = {}
    for this_h,this_ht,this_b in zip(h,ht,b):
        if this_h and this_ht:
            this_covs,this_cts = numpy.histogram([random.choice(range(this_ht)) for i in range(this_h)],bins=range(this_ht))
            coverage.update(zip(['%s.%s' % (this_b,this_ct) for this_ct in this_cts],this_covs))
    return coverage

def simulate_readcounts(num_reads,num_indiv,stdev_ratio):
    reads_mean = num_reads/num_indiv
    reads_std = reads_mean/3
    readcounts = [int(i) for i in numpy.random.normal(reads_mean,reads_std,size=num_indiv)]
    return readcounts

def simulate_indiv_coverage(num_reads,num_indiv,stdev_ratio,frags,mean,stdev,lower_bound,upper_bound):
    print >> sys.stderr, 'generate %s individuals, %s reads target ...' % (num_indiv, num_reads),
    readcounts = simulate_readcounts(num_reads,num_indiv,stdev_ratio)
    print >> sys.stderr, 'actual total generated reads: %s' % sum(readcounts)
    indiv_coverage = {}
    for i,rct in enumerate(readcounts):
        print >> sys.stderr, 'generate individual %s, %s reads' % (i,rct)
        indiv_coverage[(i,rct)] = simulate_coverage(frags,mean,stdev,rct,lower_bound,upper_bound)
    return indiv_coverage
    
def calc_shared_matrix(indiv_coverage,dp_cut):
    shared_mat = numpy.zeros((len(indiv_coverage),len(indiv_coverage)),dtype=int)
    for (i1,rct1),ind1_cov in sorted(indiv_coverage.items()):
        print >> sys.stderr, i1,rct1,dp_cut
        for locus in ind1_cov.keys():
            if ind1_cov.get(locus,0) >= dp_cut:
                for (i2,rct2),ind2_cov in sorted(indiv_coverage.items()):
                    if ind2_cov.get(locus,0) >= dp_cut:
                        shared_mat[i1,i2] += 1

    return shared_mat

def run_full_plot(frags,clust_members,indivs_by_pool,pool,samp_sd,sf1=0.95,sf2=0.9):
    from pylab import plot, scatter, subplot, savefig, clf, title
    import iplot, numpy
    matched_reads = []
    losd_reads = []
    hisd_reads = []

    for ind in indivs_by_pool[pool]:
        numreads = sum(clust_members[ind].values())
        print >> sys.stderr, ind,numreads
        mr = simulate_coverage(frags,300, samp_sd,numreads,0,1000)
        matched_reads.append(mr)
        mr = simulate_coverage(frags,300, samp_sd-4,numreads,0,1000)
        losd_reads.append(mr)
        mr = simulate_coverage(frags,300, samp_sd+4,numreads,0,1000)
        hisd_reads.append(mr)

    matched_reads_d = dict([((i,sum(mr.values())),mr) for i,mr in enumerate(matched_reads)])
    losd_reads_d = dict([((i,sum(mr.values())),mr) for i,mr in enumerate(losd_reads)])
    hisd_reads_d = dict([((i,sum(mr.values())),mr) for i,mr in enumerate(hisd_reads)])
    obs_reads = [clust_members[ind] for ind in indivs_by_pool[pool]]
    obs_reads_d = dict([((i,sum(mr.values())),mr) for i,mr in enumerate(obs_reads)]) 

    clf()
    subplot(2,2,1)
    plotcols = iplot.subspectrum(len(indivs_by_pool[pool]))
    for col,ind,mr in zip(plotcols, indivs_by_pool[pool], matched_reads):
        plot(numpy.log1p(sorted(clust_members[ind].values(),reverse=True))[10:100000],color=col)
        plot(numpy.log1p(sorted(mr.values(),reverse=True))[10:100000],color=col,ls='--',lw=1)

    title(pool)
    subplot(2,2,2)
    sf = 1
    x,y = Util.dezip(sorted([(sum(v.values()),numpy.mean([i for i in v.values() if i>1])) for v in matched_reads]))
    plot(x,sf*numpy.array(y),'k')
    x,y = Util.dezip(sorted([(sum(v.values()),numpy.mean([i for i in v.values() if i>1])) for v in losd_reads]))
    plot(x,sf*numpy.array(y),'k:')
    x,y = Util.dezip(sorted([(sum(v.values()),numpy.mean([i for i in v.values() if i>1])) for v in hisd_reads]))
    plot(x,sf*numpy.array(y),'k--')
    x,y = Util.dezip([(sum(clust_members[ind].values()),numpy.mean([i for i in clust_members[ind].values() if i>1])) for ind in indivs_by_pool[pool]])
    scatter(x,y,c='none',lw=2)

    subplot(2,2,3)
    sf = sf1
    x,y = Util.dezip([(sum(clust_members[ind].values()),len([i for i in clust_members[ind].values() if i >= 4])) for ind in indivs_by_pool[pool]])
    scatter(x,y,c='none',edgecolor='r',lw=2)
    x,y = Util.dezip(sorted([(sum(v.values()),len([i for i in v.values() if i >= 4])) for v in matched_reads]))
    plot(x,sf*numpy.array(y),'r')
    x,y = Util.dezip(sorted([(sum(v.values()),len([i for i in v.values() if i >= 4])) for v in losd_reads]))
    plot(x,sf*numpy.array(y),'r:')
    x,y = Util.dezip(sorted([(sum(v.values()),len([i for i in v.values() if i >= 4])) for v in hisd_reads]))
    plot(x,sf*numpy.array(y),'r--')

    x,y = Util.dezip([(sum(clust_members[ind].values()),len([i for i in clust_members[ind].values() if i >= 7])) for ind in indivs_by_pool[pool]])
    scatter(x,y,c='none',edgecolor='g',lw=2)
    x,y = Util.dezip(sorted([(sum(v.values()),len([i for i in v.values() if i >= 7])) for v in matched_reads]))
    plot(x,sf*numpy.array(y),'g')
    x,y = Util.dezip(sorted([(sum(v.values()),len([i for i in v.values() if i >= 7])) for v in losd_reads]))
    plot(x,sf*numpy.array(y),'g:')
    x,y = Util.dezip(sorted([(sum(v.values()),len([i for i in v.values() if i >= 7])) for v in hisd_reads]))
    plot(x,sf*numpy.array(y),'g--')

    x,y = Util.dezip([(sum(clust_members[ind].values()),len([i for i in clust_members[ind].values() if i >= 10])) for ind in indivs_by_pool[pool]])
    scatter(x,y,c='none',edgecolor='b',lw=2)
    x,y = Util.dezip(sorted([(sum(v.values()),len([i for i in v.values() if i >= 10])) for v in matched_reads]))
    plot(x,sf*numpy.array(y),'b')
    x,y = Util.dezip(sorted([(sum(v.values()),len([i for i in v.values() if i >= 10])) for v in losd_reads]))
    plot(x,sf*numpy.array(y),'b:')
    x,y = Util.dezip(sorted([(sum(v.values()),len([i for i in v.values() if i >= 10])) for v in hisd_reads]))
    plot(x,sf*numpy.array(y),'b--')

    subplot(2,2,4)
    sf = sf2
    shared_mat = calc_shared_matrix(matched_reads_d,4)
    lo_shared_mat = calc_shared_matrix(losd_reads_d,4)
    hi_shared_mat = calc_shared_matrix(hisd_reads_d,4)
    obs_shared_mat = calc_shared_matrix(obs_reads_d,4)
    x,y = Util.dezip(sorted([(k[1],numpy.mean(obs_shared_mat[k[0]])) for k in obs_reads_d.keys()]))
    scatter(x,y,c='none',edgecolor='r',lw=2)
    x,y = Util.dezip(sorted([(k[1],numpy.mean(shared_mat[k[0]])) for k in matched_reads_d.keys()]))
    plot(x,sf*numpy.array(y),'r')
    x,y = Util.dezip(sorted([(k[1],numpy.mean(lo_shared_mat[k[0]])) for k in losd_reads_d.keys()]))
    plot(x,sf*numpy.array(y),'r:')
    x,y = Util.dezip(sorted([(k[1],numpy.mean(hi_shared_mat[k[0]])) for k in hisd_reads_d.keys()]))
    plot(x,sf*numpy.array(y),'r--')

    shared_mat = calc_shared_matrix(matched_reads_d,7)
    lo_shared_mat = calc_shared_matrix(losd_reads_d,7)
    hi_shared_mat = calc_shared_matrix(hisd_reads_d,7)
    obs_shared_mat = calc_shared_matrix(obs_reads_d,7)
    x,y = Util.dezip(sorted([(k[1],numpy.mean(obs_shared_mat[k[0]])) for k in obs_reads_d.keys()]))
    scatter(x,y,c='none',edgecolor='g',lw=2)
    x,y = Util.dezip(sorted([(k[1],numpy.mean(shared_mat[k[0]])) for k in matched_reads_d.keys()]))
    plot(x,sf*numpy.array(y),'g')
    x,y = Util.dezip(sorted([(k[1],numpy.mean(lo_shared_mat[k[0]])) for k in losd_reads_d.keys()]))
    plot(x,sf*numpy.array(y),'g:')
    x,y = Util.dezip(sorted([(k[1],numpy.mean(hi_shared_mat[k[0]])) for k in hisd_reads_d.keys()]))
    plot(x,sf*numpy.array(y),'g--')

    shared_mat = calc_shared_matrix(matched_reads_d,10)
    lo_shared_mat = calc_shared_matrix(losd_reads_d,10)
    hi_shared_mat = calc_shared_matrix(hisd_reads_d,10)
    obs_shared_mat = calc_shared_matrix(obs_reads_d,10)
    x,y = Util.dezip(sorted([(k[1],numpy.mean(obs_shared_mat[k[0]])) for k in obs_reads_d.keys()]))
    scatter(x,y,c='none',edgecolor='b',lw=2)
    x,y = Util.dezip(sorted([(k[1],numpy.mean(shared_mat[k[0]])) for k in matched_reads_d.keys()]))
    plot(x,sf*numpy.array(y),'b')
    x,y = Util.dezip(sorted([(k[1],numpy.mean(lo_shared_mat[k[0]])) for k in losd_reads_d.keys()]))
    plot(x,sf*numpy.array(y),'b:')
    x,y = Util.dezip(sorted([(k[1],numpy.mean(hi_shared_mat[k[0]])) for k in hisd_reads_d.keys()]))
    plot(x,sf*numpy.array(y),'b--')


def mean_cov_over_depths(frags,sizing_mean,sizing_std,nreads,lower_bound,upper_bound,run_reps,depths):

    thisdict = defaultdict(list)
    for rp in range(run_reps):
        mr = simulate_coverage(frags,sizing_mean,sizing_std,nreads,lower_bound,upper_bound)
        for dp in depths:
            thisdict[dp].append(len([i for i in mr.values() if i >= dp]))

    return dict([(k,numpy.mean(v)) for k,v in thisdict.items()])

def append_to_dict_val(fulldict,thisdict,nreads,done):
    for k,v in thisdict.items():
        if not k in done:
            fulldict[k]['y'].append(v)
            fulldict[k]['x'].append(nreads)

def slope(p1,p2):
    return float(p2[1]-p1[1])/float(p2[0]-p1[0])


def derivative(x,y,pad=True):
    xy = zip(x,y)
    if pad:
        x0,y0 = Util.dezip([(x[i],slope(xy[i],xy[i+1])) for i in range(len(xy)-1)]+[(x[-1],slope(xy[-2],xy[-1]))])
    else:
        x0,y0 = Util.dezip([(x[i],slope(xy[i],xy[i+1])) for i in range(len(xy)-1)])
    return x0,y0

def plot_derivatives(x,y,subpl,xmax, sm=20):
    x,y = x[:-sm/2],Util.smooth(y,sm)[:-sm/2]
    x0,y0 = derivative(x,y,False)
    x0,y0 = x0[:-sm/2],Util.smooth(y0,sm)[:-sm/2]
    x1,y1 = derivative(x0,y0,False)
    x1,y1 = x1[:-sm/2],Util.smooth(y1,sm)[:-sm/2]

    pylab.subplot(*subpl)
    #pylab.cla()
    pylab.plot(x,y,'b',label='fn')
    pylab.twinx()
    pylab.plot(x0,y0,'g',label='d0')
    pylab.yticks([])
    pylab.twinx()
    pylab.plot(x1,y1,'r',label='d1')
    pylab.yticks([])
    pylab.plot([x1[0],x1[-1]],[0,0],'k:')
    pylab.xlim(0,xmax)


    

def find_read_coverage_inflection(frags,sizing_mean,sizing_std,depths=[7,10],min_reads=10000,max_reads=5000000,reads_step=10000,run_reps=10,lower_bound=0,upper_bound=1000,run_past_max=100,draw_plot=False,sm=20):
    resultdict = {}
    for dp in depths:
        resultdict[dp] = {'x':[],'y':[]}

    done = set([])

    run_past = 0
    if draw_plot:
        import pylab
        pylab.figure(1)
        print >> sys.stderr, 'plotting to %s' % draw_plot

    nreads = min_reads
    append_to_dict_val(resultdict,mean_cov_over_depths(frags,sizing_mean,sizing_std,nreads,lower_bound,upper_bound,run_reps,depths),nreads,done)
    print >> sys.stderr, 'begin iteration'
    for nreads in xrange(min_reads+reads_step,max_reads,reads_step):
        #check run_over
        if run_past_max and max(depths) in done:
            run_past += 1
        if run_past_max and run_past > run_past_max:
            break
        print >> sys.stderr, nreads
        append_to_dict_val(resultdict,mean_cov_over_depths(frags,sizing_mean,sizing_std,nreads,lower_bound,upper_bound,run_reps,depths),nreads,done)
        if len(resultdict.values()[0]['x']) > (2*sm)+3:
            for i,dp in enumerate(depths):
                if dp in done: continue
                x,y = resultdict[dp]['x'], resultdict[dp]['y']
                x,y = x[:-sm/2],Util.smooth(y,sm)[:-sm/2]
                x0,y0 = derivative(x,y,False)
                x0,y0 = x0[:-sm/2],Util.smooth(y0,sm)[:-sm/2]
                x1,y1 = derivative(x0,y0,False)
                x1,y1 = x1[:-sm/2],Util.smooth(y1,sm)[:-sm/2]

                if draw_plot:
                    pylab.subplot(len(depths),1,i+1)
                    pylab.cla()
                    pylab.plot(x,y)
                    pylab.twinx()
                    pylab.plot(x0,y0)
                    pylab.yticks([])
                    pylab.twinx()
                    pylab.plot(x1,y1)
                    pylab.plot([x1[0],x1[-1]],[0,0],'k:')
                    pylab.savefig(draw_plot)
                    
                    
                  
                #print >> sys.stderr, dp,x,y,x0,y0,x1,y1
                if y0[-1] < max(y0) and y1[-1] < 0:
                    print >> sys.stderr, 'check %sx covergence ...' % dp,
                    y1a = numpy.array(y1[numpy.argmax(y1):])
                    convg = y1a[y1a<0] ** 2
                    if convg[-1] < ( 0.01 * max(convg) ):
                        done.add(dp)
                        print >> sys.stderr, 'inflection %sx: %s reads, %s regions, store as %s %s' % (dp,x[-1],y[-1],x[-7],y[-7])
                        resultdict[dp]['x'] = resultdict[dp]['x'][:-(sm+2)]
                        resultdict[dp]['y'] = resultdict[dp]['y'][:-(sm+2)]
                    else:
                        print >> sys.stderr, 'incomplete'
        if done == set(depths):
            print >> sys.stderr, 'all converged'
            break
            
                    
        
    return resultdict


if __name__ == '__main__':

    genome, site1, site2, sizing_mean, sizing_std, cov_cuts, sm, reads_step, outroot, draw_plot = sys.argv[1:]

    sizing_mean = float(sizing_mean)
    sizing_std = float(sizing_std)
    cov_cuts = [int(i) for i in cov_cuts.split(',')]

    # a few paramters presumed constant
    # bounds of interesting size fragments
    lower_bound = 0
    upper_bound = 1000
    # individual read count mean / stdev_ratio
    stdev_ratio = 3

    reads_step=int(reads_step)
    sm=int(sm)
    min_reads = reads_step
    max_reads = 10000000

    run_reps = 5
    conv_reps = 5

    run_past_max = 20


    try:
        os.makedirs(outroot)
    except:
        pass

    frags_out = os.path.join(outroot,'%s_%s-%s_frags.list' % (os.path.basename(genome).split('.')[0],site1,site2))

    if os.path.exists(frags_out):
        print >> sys.stderr, '%s exists, loading ...' % frags_out,
        frags = eval(open(frags_out).read())
    else:
        print >> sys.stderr, 'load fasta %s ...' % genome,
        fa = Seq.Fasta(genome)
        print >> sys.stderr, 'simulate cuts'
        frags = cut_fa_obj(fa, site1, site2)
        del fa
        print >> sys.stderr, 'write %s ... ' % frags_out,
        open(frags_out,'w').write(frags.__repr__())
        
    print >> sys.stderr, 'done'
    
    ccstr = '-'.join([str(i) for i in cov_cuts])
    for i in range(conv_reps):
        rd_out = os.path.join(outroot,'%s_%s-%s_reads%s-%sstep%s_cov%s_sm%s_resultdict%s.dict' % (os.path.basename(genome).split('.')[0],site1,site2,min_reads,max_reads,reads_step,ccstr,sm,i))
        if eval(draw_plot):
            this_plot = os.path.join(outroot,'%s_%s-%s_reads%s-%sstep%s_cov%s_sm%s_resultdict%s.png' % (os.path.basename(genome).split('.')[0],site1,site2,min_reads,max_reads,reads_step,ccstr,sm,i))
        else:
            this_plot = False
        if not os.path.exists(rd_out):
            #uncomment for "leapfrog" parallelization:
            os.system('touch '+rd_out)
            rd = find_read_coverage_inflection(frags,sizing_mean,sizing_std,cov_cuts,min_reads=min_reads,reads_step=reads_step,max_reads=max_reads,run_reps=run_reps,run_past_max=run_past_max,sm=sm,draw_plot=this_plot)
            print >> sys.stderr, 'write resutdict to %s' % rd_out
            open(rd_out,'w').write(rd.__repr__())

