#!/usr/bin/env python
'''
fit_gaussian_to_cluster_rep.py mean stdev_lo stdev_hi cluster_rep source_fragments sample_dists sample_dist_draws log
'''

import numpy, random, os, sys

def sample_from_fragments(sample_hist,fragments_hist):
    '''
    draws a sample representation distribution for a given fragment distribution and sample target distribution, i.e. gaussian
    '''
    coverage = []

    for this_hs,this_hf in zip(sample_hist,fragments_hist):
        if this_hs and this_hf:
            coverage.extend((numpy.histogram([random.choice(range(this_hf)) for i in range(this_hs)],bins=range(this_hf))[0]))

    return sorted(coverage,reverse=True)

if __name__ == '__main__':

    sizemin = 0
    sizemax = 1000

    mean,stdl,stdh,stdstep,clustrepf,sourcefragf,sdists,sdistdraws,lfn,outfile = sys.argv[1:]
    
    mean,stdl,stdh,stdstep,sdists,sdistdraws = [float(i) for i in [mean,stdl,stdh,stdstep,sdists,sdistdraws]]

    logfn = eval(lfn)
    if logfn == None: logfn = lambda x: x

    print >> sys.stderr, 'load source fragments from %s ...' % sourcefragf, 
    sourcefrags = eval(open(sourcefragf).read())
    print >> sys.stderr, '%s fragments, total length %s' % (len(sourcefrags), sum(sourcefrags))

    hf,bf = numpy.histogram(sourcefrags,bins=range(sizemin,sizemax))
    
    print >> sys.stderr, 'load cluster representation from %s ...' % clustrepf,
    clustreps = eval(open(clustrepf).read())
    print >> sys.stderr, '%s clusters, total reads %s' % (len(clustreps), sum(clustreps))
    
    result = {}

    for stdv in numpy.arange(stdl,stdh,stdstep):
        result[stdv] = {}
        for sdistnum in range(int(sdists)):
            print >> sys.stderr, 'stdev %s draw %s' % (stdv,sdistnum),
            result[stdv][sdistnum] = []
            sample = numpy.random.normal(mean,stdv,size=sum(clustreps))
            hs,bs = numpy.histogram(sample,bins=range(0,1000))
            for sdistdrawnum in range(int(sdistdraws)):
                scov = sample_from_fragments(hs,hf)
                try:
                    rsq = numpy.corrcoef(logfn(clustreps[:len(scov)]),logfn(scov))[0][1]
                    result[stdv][sdistnum].append(rsq)
                    print >> sys.stderr, '.',
                except:
                    rsq = numpy.corrcoef(logfn(clustreps),logfn(scov[:len(clustreps)]))[0][1]
                    result[stdv][sdistnum].append(rsq)
                    print >> sys.stderr, 'x',
                                    
                    
            print >> sys.stderr,numpy.median(result[stdv][sdistnum])

    
    open(outfile,'w').write(result.__repr__())
    
    
