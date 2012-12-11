import Util,sys
from short_read_analysis import variant_detection
from numpy import prod
from collections import defaultdict

def init_data(l3_li):
    minor_hom = [] #commited to minor allele homozygote
    minor_het = [] #commited to minor allele heterozygote
    hetsort = sorted(l3_li,key = lambda x:x[0][1]) #sort by best het 
    af = {}
    #L maf = 0
    af[(len(minor_het)) + (len(minor_hom)*2)] = prod([el[0][0] for el in hetsort]) * prod([el[0][1] for el in minor_het]) *  prod([el[0][2] for el in minor_hom])
    #start with the most likely het
    minor_het.append(hetsort.pop())
    return hetsort,minor_het,minor_hom,af


def af_method1(hetsort,minor_het,minor_hom,af,verbose=True):
    # METHOD 1:
    # allows all possible category reassignments in calculation; I think this is correct

    while hetsort or minor_het: #while the "homozygous reference" category still has genotypes
        #record the current breakdown
        af[(len(minor_het)) + (len(minor_hom)*2)] = prod([el[0][0] for el in hetsort]) * prod([el[0][1] for el in minor_het]) *  prod([el[0][2] for el in minor_hom]) 
        if hetsort:
            homsort = sorted(hetsort,key = lambda x:x[0][2]) #copy sorted by best homozygote
            hom2hom = prod([el[0][0] for el in homsort[:-1]]) \
            * prod([el[0][1] for el in minor_het[:-1]]) \
            * prod([el[0][2] for el in minor_hom]) \
            * (minor_het and minor_het[-1][0][0] or 0) \
            * homsort[-1][0][2]
            hom2het = prod([el[0][0] for el in hetsort[:-1]]) \
            * prod([el[0][1] for el in minor_het]) \
            * prod([el[0][2] for el in minor_hom]) \
            * hetsort[-1][0][1]
        else:
            hom2hom = 0
            hom2het = 0
        if minor_het:
            homsort_minor_het = sorted(minor_het,key = lambda x:x[0][2])
            het2hom = prod([el[0][0] for el in hetsort]) \
            * prod([el[0][1] for el in homsort_minor_het[:-1]]) \
            * prod([el[0][2] for el in minor_hom]) \
            * homsort_minor_het[-1][0][2]
        else:
            het2hom = 0
        if hom2hom > max(hom2het,het2hom) and minor_het:
            if verbose: print >>sys.stderr,'%s beats %s, %s move %s to hom nonref, move %s to hom ref' % \
            (hom2hom,hom2het,het2hom,homsort[-1],minor_het[-1])
            if homsort:
                minor_hom.append(homsort.pop())
            if minor_het:
                homsort.append(minor_het.pop())
            hetsort = sorted(homsort,key = lambda x:x[0][1])
        #otherwise check if take the next het, as above
        elif hom2het > max(hom2hom,het2hom):
            if verbose: print >>sys.stderr,'%s beats %s, %s move %s to het' % (hom2het,hom2hom,het2hom,hetsort[-1])
            minor_het.append(hetsort.pop())
        #otherwise "promote" a het
        elif het2hom > max(hom2hom,hom2het):
            if verbose: print >>sys.stderr,'%s beats %s, %s move %s to hom nonref' % (het2hom,hom2hom,hom2het,homsort_minor_het[-1])
            minor_hom.append(homsort_minor_het.pop())
            minor_het = sorted(homsort_minor_het,key = lambda x:x[0][1])
        else:
            if not hetsort: #no hom ref, promote a het
                if verbose: print >>sys.stderr,'no winner, no hom ref, move %s to hom nonref' % (homsort_minor_het[-1].__repr__())
                minor_hom.append(homsort_minor_het.pop())
                minor_het = sorted(homsort_minor_het,key = lambda x:x[0][1])
            elif not minor_het: #no het, promote a hom ref
                if verbose: print >>sys.stderr,'no winner, no minor_het, move %s to het' % (hetsort[-1].__repr__())
                minor_het.append(hetsort.pop())
            elif homsort_minor_het[-1][0][2] > hetsort[-1][0][1]: # het -> hom nonref beats hom ref -> het
                if verbose: print >>sys.stderr,'no winner, het to hom nonref beats hom ref to het move %s to hom nonref' % (homsort_minor_het[-1].__repr__())
                minor_hom.append(homsort_minor_het.pop())
                minor_het = sorted(homsort_minor_het,key = lambda x:x[0][1])
            else:
                if verbose: print >>sys.stderr,'no winner, hom ref to het beats het to hom nonref, move %s to het' % (hetsort[-1].__repr__())
                minor_het.append(hetsort.pop())
    return af

def af_method2(hetsort,minor_het,minor_hom,af):
    ############################################################
    # METHOD 2
    # considers only A) hom ref to het or B) hom ref to hom nonref 
    # (with accompanying het to hom ref to balance)
    #
    # guaranteed (I think) to consider all allele frequencies, 
    # but can get stuck in local minima if het to hom nonref is higher likelihood than A or B above


    while hetsort:
        af[(len(minor_het)) + (len(minor_hom)*2)] = prod([el[0][0] for el in hetsort]) * prod([el[0][1] for el in minor_het]) *  prod([el[0][2] for el in minor_hom])
        homsort = sorted(hetsort,key = lambda x:x[0][2])
        if prod([el[0][0] for el in hetsort]) * prod([el[0][1] for el in minor_het[:-1]]) *  prod([el[0][2] for el in minor_hom]) * minor_het[-1][0][0] * homsort[-1][0][2] > prod([el[0][0] for el in hetsort]) * prod([el[0][1] for el in minor_het]) *  prod([el[0][2] for el in minor_hom]) * hetsort[-1][0][1]:
            minor_hom.append(homsort.pop())
            homsort.append(minor_het.pop())
            hetsort = sorted(homsort,key = lambda x:x[0][1])
        else:
            minor_het.append(hetsort.pop())
            

    while minor_het:
        af[(len(minor_het)) + (len(minor_hom)*2)] = prod([el[0][0] for el in hetsort]) * prod([el[0][1] for el in minor_het]) *  prod([el[0][2] for el in minor_hom])
        minor_hom.append(minor_het.pop())

    return af

def cut_fn(sd): #load sites with at least one non-ref site, 10 genotyped individuals
    return len(sd['indiv_gt']) > 10 and sd['mac'] >= 1

if __name__ == "__main__":

    #use the above, plus toss any calls with max quality < 4 (GQ)
    vcf = variant_detection.load_vcf('/n/hoekstrafs1/test-stampy/110910-lane5_stampy.vcf', cutoff_fn=cut_fn, indiv_gt_phred_cut=4)

    dephred = lambda x: 10**(x/-10.)

    l3_li = []
    for v in vcf.values()[9]['indiv_gt'].values():
        l3 = [dephred(int(p)) for p in v['PL'].split(',')]
        l3_li.append(([l/sum(l3) for l in l3],v))
       

    hetsort,minor_het,minor_hom,af = init_data(l3_li)

    af = af_method1(hetsort,minor_het,minor_hom,af)
           
    plot(*Util.dezip(sorted(af.items())))
    len(af) == 2*len(l3_li) #can be false, but so far all missing L = 0.0 (so irrelevant)


    hetsort,minor_het,minor_hom,af2 = init_data(l3_li)       

    af2 = af_method2(hetsort,minor_het,minor_hom,af2)

    plot(*Util.dezip(sorted(af2.items())))
    len(af2) == 2*len(l3_li) #always true so far (I think this is right)

    af_tot = defaultdict(float)
    for i in range(1000):
        l3_li = []
        for v in vcf.values()[i]['indiv_gt'].values():
            l3 = [dephred(int(p)) for p in v['PL'].split(',')]
            l3_li.append(([l/sum(l3) for l in l3],v))
        if len(l3_li[0][0]) != 3: continue
        hetsort,minor_het,minor_hom,af = MLaf.init_data(l3_li)
        af = MLaf.af_method1(hetsort,minor_het,minor_hom,af,False)
        af_expand = dict([(int(k* ((num_ind*2)/float(len(af)))),v) for k,v in af.items()])
        af_exp_fill = dict(zip(range((num_ind*2)),Util.interpolate_missing_values([af_expand.get(i,None) for i in range((num_ind*2))])))
        for k,v in af_exp_fill.items():
            if v > 0:
                af_tot[k]+=log10(v)

    plot(*Util.dezip(sorted(af_tot.items())[:-1])) #last element is broken; need to figure out edge effects...



    plot(*Util.dezip(sorted(af_tot.items())))



