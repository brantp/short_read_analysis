from collections import defaultdict
import sys, numpy
import os, sys, re, Util
import Seq
import gzip

from subprocess import Popen,PIPE

def smartopen(filename,*args,**kwargs):
    '''opens with open unless file ends in .gz, then use gzip.open

    in theory should transparently allow reading of files regardless of compression
    '''
    if filename.endswith('.gz'):
        return gzip.open(filename,*args,**kwargs)
    else:
        return open(filename,*args,**kwargs)



def load_vcf(vcf,cutoff_fn=None,ding_on=100000,store_only=None,indiv_gt_phred_cut=None,store_indiv_prefix=None,drop_indiv=None,keep_indiv=None,write_thresholded_vcf=None,write_fasta_base=None,ref_in_fasta=False,regions=None,multiallelic_sites='skip',noGQ='skip'):
    '''populates and returns a site:properties dict from vcf file

    if <store_only> is set, must be a list of fields to retain

    if <cutoff_fn> is set, each sd (dict of fully parsed line, subject to store_only filter) is passed to cutoff_fn
    must be callable, and truth value of the return determines retention of that site

    <multiallelic_sites> : ["skip"|"keep"|FUNCTION] (if FUNCTION, will be called with sd and must return either a fixed sd or None)
    '''

    ma_sites_processed = 0
    ma_sites_kept = 0

    if noGQ == 'binom':
        print >> sys.stderr, 'noQG is binom; this is a very sad hack based on the observation:'
        print >> sys.stderr, 'phred( P(no variant reads|variant>20% in library) ) == number of reads'
        print >> sys.stderr, 'PLEASE MAKE SURE YOU KNOW WHAT YOU ARE DOING'

    if write_thresholded_vcf is not None:
        ofh = smartopen(write_thresholded_vcf,'w')

    if write_fasta_base is not None:
        import Seq
        xd = { '0':'REF','1':'ALT' }
        faroot = os.path.dirname(write_fasta_base)
        try:
            os.makedirs(faroot)
        except:
            pass
        faseqs = {}
        this_chrom = None
        last_data = 0
        indivs = []

    vcf_data = {}
    
    i = 0
    for line in smartopen(vcf):
        if write_thresholded_vcf is not None and line.startswith('##'):
            ofh.write(line)
        
        if i % ding_on == 0: print >> sys.stderr, 'reading',i
        i += 1

        #skip bogus lines -- probably want to warn on this somehow...
        if '\x00' in line:
            print >> sys.stderr, 'Null byte in line %s, skipping' % i
            continue

        if line.startswith('#CHROM'):
            if write_thresholded_vcf is not None:
                ofh.write(line)
            headers = line[1:].split()
            exp_elements = len(line.split())
            FORMAT = headers.index('FORMAT')
            indivs = headers[FORMAT+1:]
            if write_fasta_base is not None: #DON'T INCLUDE DROP_INDIV
                if drop_indiv is not None:
                    indivs = list(set(headers[FORMAT+1:]) - set(drop_indiv))
                else:
                    indivs = headers[FORMAT+1:]
            if drop_indiv is not None or keep_indiv is not None:
                all_keep = set(indivs).intersection(set(keep_indiv is None and indivs or keep_indiv)) - set(drop_indiv and drop_indiv or [])
                all_drop = set(indivs) - all_keep
                if len(all_keep) <= len(all_drop):
                    #test keeps
                    print >> sys.stderr, 'length of keep set (%s) less than length of drop set (%s); use keep' % (len(all_keep),len(all_drop))
                    keep_indiv = list(all_keep)
                    drop_indiv = None
                else:
                    #test drops
                    print >> sys.stderr, 'length of keep set (%s) greater than length of drop set (%s); use drop' % (len(all_keep),len(all_drop))
                    keep_indiv = None
                    drop_indiv = list(all_drop)
        elif line.startswith('#'):
            continue
        else:
            #extract site stats
            fields = line.split()
            if len(fields) != exp_elements:
                print >>sys.stderr, 'unexpected length, line %s (exp %s obs %s)' % (i,exp_elements,len(fields))
                continue

            if regions is not None:
                if not fields[0] in regions:
                    continue

            if write_fasta_base is not None:
                if this_chrom != fields[0]: #moving on, write last
                    if this_chrom is not None and last_data != 0:
                        faofh = open(write_fasta_base+this_chrom+'.fa','w')
                        for ind,slist in faseqs.items():
                            if not set(slist).issubset(set(['N','n','R', 'S', 'B', 'D', 'W', 'V', 'Y', 'H', 'K', 'M'])):
                                faofh.write('>%s\n%s\n' % (ind, ''.join(slist)))
                        faofh.close()
                    this_chrom = fields[0]
                    faseqs = {}
                    for ind in indivs: 
                        faseqs[ind] = []
                    if ref_in_fasta:
                        faseqs['REF'] = []
                    last_data = 0
            
            sd = dict(zip(headers[:FORMAT],fields[:FORMAT]))
            key = (sd['CHROM'],sd['POS'])
            try:
                infostr = sd.pop('INFO')
                sd.update(dict([el.split('=') for el in infostr.split(';') if '=' in el]))
            except KeyError:
                pass

            sd['indiv_gt'] = {}
            for ind,gt in zip(headers[FORMAT+1:],fields[FORMAT+1:]):
                if store_indiv_prefix is not None and not ind.startswith(store_indiv_prefix): continue
                #drop_indiv HERE
                if drop_indiv is not None and ind in drop_indiv: continue
                if keep_indiv is not None and not ind in keep_indiv: continue
                if ':' in gt:
                    this_gt = dict(zip(fields[FORMAT].split(':'),gt.split(':')))
                    if this_gt.has_key('GQ'):
                        if indiv_gt_phred_cut is None or float(this_gt['GQ']) >= indiv_gt_phred_cut:
                            sd['indiv_gt'][ind] = this_gt
                    else:
                        if noGQ == 'keep':
                            sd['indiv_gt'][ind] = this_gt #store anyway
                        elif noGQ == 'skip':
                            continue #move on to next indiv
                        elif noGQ == 'binom': #MEGAHACK -- phred( P(no variant reads|variant>20% in library) ) == number of reads
                            this_gt['GQ'] == this_gt['DP']
                        else:
                            print >> sys.stderr, 'GQ field absent (%s) and noGQ parameter not in [keep|skip] (is: %s)' % (gt,noGQ)
                            raise

            if write_thresholded_vcf is not None: #ugly hack to "zero out" below-thresh GTs
                for i in range(FORMAT+1,len(fields)):
                    if ':' in fields[i]:
                        this_gt = dict(zip(fields[FORMAT].split(':'),fields[i].split(':')))
                        if indiv_gt_phred_cut is not None and float(this_gt['GQ']) < indiv_gt_phred_cut:
                            fields[i] = './.'

            if len(sd['indiv_gt']) > 0:
                if ',' in sd['ALT']:
                    ma_sites_processed += 1
                    if multiallelic_sites == 'skip':
                        continue
                    elif multiallelic_sites == 'keep':
                        ma_sites_kept += 1
                    else:
                        try:
                            sd = multiallelic_sites(sd)
                        except:
                            print >> sys.stderr, 'multiallelic sites resolution call failed at %s,%s' % (sd['CHROM'],sd['POS']) 
                            sd = None
                        if sd is None:
                            continue
                        ma_sites_kept += 1
                else:
                    sd = add_custom_fields(sd)
            else:
                continue #this is questionable -- ever want to proceed with a genotype-less site?
                
            #before we toss data, build sequence if desired
            if write_fasta_base is not None:
                pad_n = int(sd['POS']) - last_data - 1
                if cutoff_fn is None or cutoff_fn(sd):
                    var_accept = True
                else:
                    var_accept = False

                if ref_in_fasta:
                    faseqs['REF'].extend(['N']*pad_n)
                    faseqs['REF'].append(sd['REF'])
                
                for ind in indivs:
                    faseqs[ind].extend(['N']*pad_n)
                    if sd['indiv_gt'].has_key(ind):
                        gt_li = sd['indiv_gt'][ind]['GT'].split('/')
                        #if '1' not in gt_li or var_accept:
                        if var_accept: #only retain data if whole column is valid; contrast with retaining invariants even in "reject" sites ('1' not in gt_li)
                            this_gt = Seq.degenerate_iupac([[sd[xd[allele]] for allele in gt_li]]) #FIX!
                        else:
                            this_gt = 'n'
                    else:
                        this_gt = 'N'
                    faseqs[ind].append(this_gt)
                last_data = int(sd['POS'])
                                           

            if cutoff_fn is None or cutoff_fn(sd):
                if store_only is not None:
                    keep_sd = {}
                    for k in sd:
                        if k in store_only:
                            keep_sd[k] = sd[k]
                    sd = keep_sd
                if len(sd) > 0:
                    vcf_data[key] = sd
                if write_thresholded_vcf is not None:
                    ofh.write('\t'.join(fields) + '\n')

                    
    if write_thresholded_vcf is not None:           
        ofh.close()

    if write_fasta_base is not None:
        faofh = open(write_fasta_base+this_chrom+'.fa','w')
        for ind,slist in faseqs.items():
            faofh.write('>%s\n%s\n' % (ind, ''.join(slist)))
        faofh.close()

    print >> sys.stderr, '%s multi-allelic sites processed; %s retained' % (ma_sites_processed,ma_sites_kept)

    return vcf_data

def add_custom_fields(sd):
    sd['hwe'] = hwe(sd['indiv_gt'])
    sd['mac'] = maf(sd['indiv_gt'])
    sd['maf'] = sd['mac'] / (2. * len(sd['indiv_gt']))
    if sd['hwe'] is None:
        sd['mafbyhwe'] = None
    else:
        sd['mafbyhwe'] = sd['maf'] / numpy.log1p(sd['hwe'])
    sd['fh'] = fract_het(sd['indiv_gt'])
    sd['totcov'] = sum([int(gtd['DP']) for gtd in sd['indiv_gt'].values()])
    sd['numind'] = len(sd['indiv_gt'])
    sd['aac'] = aac(sd['indiv_gt'])
    sd['aaf'] = sd['aac'] / (2. * len(sd['indiv_gt']))
    return sd

### functions for wrangling multiallelic calls
allelesort = lambda j,k: (int(k)*(int(k)+1)/2)+int(j)
phred = lambda x: -10 * numpy.log10(x)
dephred = lambda x: 10**(x/-10.0)

sum_allele_counts = lambda sd: reduce(lambda x,y : x+y, [numpy.array(d['AD'].split(','),dtype=int) for d in sd['indiv_gt'].values()])

def sum_invalid_allele_counts(sd):
    sac = sum_allele_counts(sd)
    sum_invalid = numpy.zeros(sac.shape,dtype=int)
    for d in sd['indiv_gt'].values():
        act = numpy.array(d['AD'].split(','),dtype=int)
        for i in list(set(d['GT'].split('/'))):
            act[int(i)] = 0
        sum_invalid += act    
    return sum_invalid

def genotype_allele_counts(sd):
    al = dict.fromkeys(range(len(sd['ALT'].split(','))+1),0)
    for d in sd['indiv_gt'].values():
        for a in d['GT'].split('/'):
            al[int(a)] += 1
    return [v for k,v in sorted(al.items())]


def reads_dropped_by_allele(sd):
    reads_dropped = dict.fromkeys(range(len(sd['ALT'].split(','))+1),0)
    for d in sd['indiv_gt'].values():
        counts = numpy.array(d['AD'].split(','),dtype=int)
        for k in list(set(d['GT'].split('/'))):
            reads_dropped[int(k)] += counts[int(k)]
    return reads_dropped

def alleleorder(nalleles):
    return sorted(list(set([tuple(sorted((str(j),str(k)))) for j in range(nalleles) for k in range(nalleles)])),key=lambda t: allelesort(*t))

def get_keep_alleles(sd):
    keep_alleles = set([str(k) for k,v in sorted(reads_dropped_by_allele(sd).items(),key= lambda x:x[1],reverse=True)[:2]])
    return keep_alleles

def get_num_alleles(sd):
    nalleles = len(sd['ALT'].split(','))+1
    return nalleles

def get_PL_dict(gtd,nalleles):
    ao = alleleorder(nalleles)
    return dict(zip(ao,map(dephred,map(float,gtd['PL'].split(',')))))

def get_keep_PL_dict(pld,keep_alleles):
    return dict([(k,v) for k,v in pld.items() if set(k).issubset(keep_alleles)])

def get_best_gt(pld):
    return max([(v,k) for k,v in pld.items()])[1]

def calc_GQ(pld,GT,GQ_method='nextbest'):
    '''methods: "sum" , "nextbest"
    '''
    if GQ_method == 'sum':
        return min(99.0,phred(1-pld[GT]/sum(pld.values())))
    elif GQ_method == 'nextbest':
        L = phred(pld[GT])
        nextbest = numpy.inf
        for this_L in pld.values():
            if nextbest > phred(this_L) > L:
                nextbest = phred(this_L)
        return min(99.0, nextbest-L)

def get_keep_gtd(gtd,nalleles,keep_alleles,old2new,new2old,GQ_method='nextbest'):
    keep_ao = alleleorder(len(keep_alleles))
    pld = get_PL_dict(gtd,nalleles)
    pld_keep = get_keep_PL_dict(pld,keep_alleles)
    gt = get_best_gt(pld_keep)
    gq = calc_GQ(pld_keep,gt,GQ_method)
    keep_ad = numpy.array(map(int,gtd['AD'].split(',')))[numpy.array(sorted(keep_alleles),dtype=int)]
    keep_gtd = {'GQ':str(gq), \
                'GT':'/'.join([old2new[a] for a in gt]), \
                'PL':','.join(map(str,[phred(pld_keep[tuple([new2old[k] for k in k_ao])]) for k_ao in keep_ao])), \
                'AD':','.join(map(str,keep_ad)), \
                'DP':str(sum(keep_ad))}                            
    return keep_gtd

def get_allele_lookup(sd):
    nalleles = get_num_alleles(sd)
    return dict(zip(map(str,range(nalleles)), [sd['REF']]+sd['ALT'].split(',')))

def biallelic_sd(sd,GQ_method='nextbest',transfer_fields=['CHROM','POS','AN','QUAL','QD','MQ0'],skip_fields=['indiv_gt']):
    new_sd = {}
    new_sd['indiv_gt'] = {}
    nalleles = get_num_alleles(sd)
    ao = alleleorder(nalleles)
    keep_alleles = get_keep_alleles(sd)

    old_alleles = get_allele_lookup(sd)
    new2old = dict(zip(map(str,range(len(keep_alleles))),sorted(keep_alleles)))
    old2new = dict(zip(sorted(keep_alleles),map(str,range(len(keep_alleles)))))
    
    delta_L = 0
    for ind,gtd in sd['indiv_gt'].items():
        new_gtd = get_keep_gtd(gtd,nalleles,keep_alleles,old2new,new2old,GQ_method)
        change = float(new_gtd['GQ']) - float(gtd['GQ'])
        delta_L += change
        new_sd['indiv_gt'][ind] = new_gtd

    ma = sum(sorted(reads_dropped_by_allele(sd).values(),reverse=True)[2:])/float(sum(sum_allele_counts(sd)))
    inv = float(sum(sum_invalid_allele_counts(sd)))/sum(sum_allele_counts(sd))
    orig_L = sum([float(gtd['GQ']) for d in sd['indiv_gt'].values()])

    new_sd['diploid_invalid_reads_fraction'] = inv
    new_sd['biallelic_invalid_reads_fraction'] = ma
    new_sd['delta_L'] = delta_L
    new_sd['orig_L'] = orig_L

    new_sd['new2old'] = new2old
    new_sd['old2new'] = old2new

    new_sd['REF'] = old_alleles[new2old['0']]
    new_sd['ALT'] = old_alleles[new2old['1']]
    new_sd['AC'] = dict(zip(map(str,range(1,nalleles)), sd['AC'].split(',')))[new2old['1']]
    new_sd['AF'] = dict(zip(map(str,range(1,nalleles)), sd['AF'].split(',')))[new2old['1']]

    for key in sd:
        if key in transfer_fields:
            new_sd[key] = sd[key]
        elif key in skip_fields:
            continue
        else:
            new_sd['orig_%s' % key] = sd[key]

    new_sd = add_custom_fields(new_sd)

    return new_sd

def resolve_multiallelic_sd_fn(min_delta_L_bi=None,max_fmulti_finv_bi=None,max_fract_invalid_reads=None,**kwargs):
    return lambda sd: resolve_multiallelic_sd(sd,min_delta_L_bi,max_fmulti_finv_bi,max_fract_invalid_reads,**kwargs)
    
def resolve_multiallelic_sd(sd,min_delta_L_bi=None,max_fmulti_finv_bi=None,max_fract_invalid_reads=None,**kwargs):
    min_err = 0.0001
    new_sd = biallelic_sd(sd,**kwargs)
    delta_L = float(new_sd['delta_L'])/float(new_sd['orig_L'])
    fract_invalid_reads = max(float(new_sd['diploid_invalid_reads_fraction']),min_err)
    fmulti_finv = float(new_sd['biallelic_invalid_reads_fraction'])/fract_invalid_reads
    if (min_delta_L_bi is not None and delta_L > min_delta_L_bi) \
           and (max_fmulti_finv_bi is not None and fmulti_finv < max_fmulti_finv_bi) \
           and (max_fract_invalid_reads is not None and fract_invalid_reads < max_fract_invalid_reads):
        return new_sd
    else:
        return None
    
    
### end multiallelic site code
        
def calc_het_llr(indiv_gtd):
    '''given a dictionary of individual genotype data per vcf_data[key]['indiv_gt'] from load_vcf

    calculates the chi-square statistic for the likelihood ratio of each site being heterozygous.

    returns a version of this dict with 'het_llr':<chi>'''

    for ind,gtd in indiv_gtd.items():
        Ls = gtd['GL']

    #nevermind, use phred-scaled P from GATK output

def hwe(indiv_gt):
    '''given an indiv_gt dict as in vcf_data

    returns the observed and expected classes for hardy-weinberg equilibrium'''

    alleles = set(reduce(lambda x,y:x+y,[re.split('[/|]',v['GT']) for v in indiv_gt.values()]))
    if len(alleles) == 1:
        return 0
    if len(alleles) > 2:
        return None

    Aid, Bid = sorted(list(alleles))

    A = 0
    a = 0

    AA = 0
    Aa = 0
    aa = 0

    for v in indiv_gt.values():
        if set(re.split('[/|]',v['GT'])) == set([Aid]):
            AA += 1
            A += 2
        elif set(re.split('[/|]',v['GT'])) == set([Aid,Bid]):
            Aa += 1
            A += 1
            a += 1
        elif set(re.split('[/|]',v['GT'])) == set([Bid]):
            aa += 1
            a += 2
    try:
        p = A/float(A+a)
    except ZeroDivisionError:
        print >> sys.stderr, indiv_gt
    q = 1-p

    tot = AA+Aa+aa

    expAA = p**2 * tot
    expAa = 2*p*q*tot
    expaa = q**2 *tot

    X2 = sum([(o-e)**2/e for o,e in zip([AA,Aa,aa],[expAA,expAa,expaa]) if e != 0])

    return X2
        
def maf(indiv_gt):
    '''returns the counts of the minor allele for an indiv_gt dict as in vcf_data above'''

    A = 0
    a = 0

    for v in indiv_gt.values():
        A += v['GT'].count('0')
        a += v['GT'].count('1')

    return min(A,a)

def aac(indiv_gt):
    '''returns the counts of the minor allele for an indiv_gt dict as in vcf_data above'''

    A = 0
    a = 0

    for v in indiv_gt.values():
        A += v['GT'].count('0')
        a += v['GT'].count('1')

    return a

def fract_het(indiv_gt):
    gts = [v['GT'] for v in indiv_gt.values()]

    return (gts.count('0/1') + gts.count('0|1') + gts.count('1|0'))/float(len(gts))

def add_vcf_data(ref_vcf,*vcfs):
    '''
    given a series of vcf files treating the first as the master, adds indiv_gt dicts via .update()
    to the relevant site dicts of the master (i.e. first listed) vcf object
    '''
    from copy import deepcopy
    ref_vcf_copy = deepcopy(ref_vcf)
    tot = len(ref_vcf_copy.keys())
    tickon = tot/20
    for i,k in enumerate(ref_vcf_copy.keys()):
        if i%tickon==0: print >> sys.stderr, '%s / %s' % (i,tot)
        for vcf in vcfs:
            if vcf.has_key(k):
                ref_vcf_copy[k]['indiv_gt'].update(vcf[k]['indiv_gt'])

    return ref_vcf_copy
        

def one_site_per_chrom(vcfdict,max_fn,filter_fn):
    '''given a vcf dictionary per load_vcf

    takes one site per chromosome (key[0]) for which filter_fn[v] == False and
    which maximizes max_fn[v]

    returns dictionary {chrom:site}
    '''

    topsites = {}
    for k,v in vcfdict.items():
        if filter_fn(v): continue #skip records that trigger filter_fn
        if topsites.has_key(k[0]):
            if max_fn(v) > max_fn(vcfdict[(k[0],topsites[k[0]])]):
                topsites[k[0]] = k[1]
        else:
            topsites[k[0]] = k[1]

    return topsites

def align_sanger_exp(seqs,reference,exp_coord_fn=None,max_dist_to_exp=3000,ind_lookup_fn=lambda x:x):
    '''
    given a fasta of sanger sequences and an alignment target, returns aligned regions as list of tuples:
    [individual,reference_start,reference_len,indiv_sequence]
    '''

    sanger_data = []
    skip = False
    rline = False
    dist_to_exp = []
    too_far = []

    for l in Popen('lastz %s[multiple] %s --ambiguous=iupac --format=maf' % (reference,seqs), shell=True, stdout=PIPE,stderr=PIPE).stdout:
        if skip == True:
            skip = False
            continue
        if l.startswith('s'):
            if rline == False:
                ref = l.split()[1]
                refstart = int(l.split()[2])
                reflen = int(l.split()[3])
                refseq = l.strip().split()[-1]
                if l.split()[4] == '-':
                    print 'skipping',l
                    skip = True
                    continue
                rline = True
            else:                               
                ind = l.split()[1].split('_')[0]
                sangSeq = l.strip().split()[-1]
                exp_pos_min = refstart - (max_dist_to_exp)
                exp_pos_max = refstart + (max_dist_to_exp)
                #print ind,refstart,reflen
                #print refseq
                #print sangSeq
                m = []                 
                for r,s in zip(refseq,sangSeq):
                    if r == s:                                                                                                                                 
                        m.append('.')                                                  
                    else:                                                                    
                        m.append(s)
                #print ''.join(m)
                if exp_coord_fn is None:
                    sanger_data.append((ind_lookup_fn(ind),ref,refstart,reflen,refseq,sangSeq,''.join(m)))
                else:
                    dist_to_exp.append((refstart - exp_coord_fn(ind)))
                    if exp_pos_min < exp_coord_fn(ind) < exp_pos_max:
                        sanger_data.append((ind_lookup_fn(ind),ref,refstart,reflen,refseq,sangSeq,''.join(m)))
                    else:
                        too_far.append((ind,refstart,exp_coord_fn(ind),len(sangSeq)))
                        
                    
                rline = False

    return sanger_data,dist_to_exp,too_far


def align_sanger(seqs,reference,exp_coord_cut_fn=None,ind_lookup_fn=lambda x:x):
    '''
    given a fasta of sanger sequences and an alignment target, returns aligned regions as list of tuples:
    [individual,reference_start,reference_len,indiv_sequence]
    '''

    sanger_data = []
    skip = False
    rline = False
    dist_to_exp = []

    for l in Popen('lastz %s[multiple] %s --ambiguous=iupac --format=maf' % (reference,seqs), shell=True, stdout=PIPE,stderr=PIPE).stdout:
        if skip == True:
            skip = False
            continue
        if l.startswith('s'):
            if rline == False:
                ref = l.split()[1]
                refstart = int(l.split()[2])
                reflen = int(l.split()[3])
                refseq = l.strip().split()[-1]
                if l.split()[4] == '-':
                    print 'skipping',l
                    skip = True
                    continue
                rline = True
            else:                               
                ind = l.split()[1].split('_')[0]
                sangSeq = l.strip().split()[-1]
                #print ind,refstart,reflen
                #print refseq
                #print sangSeq
                m = []                 
                for r,s in zip(refseq,sangSeq):
                    if r == s:
                        m.append('.')                                                  
                    else:                                                                    
                        m.append(s)
                #print ''.join(m)
                if exp_coord_cut_fn is None:
                    sanger_data.append((ind_lookup_fn(ind),ref,refstart,reflen,refseq,sangSeq,''.join(m)))
                elif exp_coord_cut_fn(ind,refstart):
                    sanger_data.append((ind_lookup_fn(ind),ref,refstart,reflen,refseq,sangSeq,''.join(m)))
                rline = False

    return sanger_data

def variants_from_sanger(sanger_data,*args):
    var_from_sanger = defaultdict(dict)
    for ind,target,refstart,reflen,refseq,sangseq,m in sanger_data:
        pos = 0
        for rbase,sbase in zip(refseq,sangseq):
            if rbase != '-':
                pos += 1
            if rbase.upper() != sbase.upper() and rbase.upper() != 'N' and sbase.upper() != 'N':
                k = (target,str(refstart+pos))
                obs_var = Seq.undegenerate_iupac(sbase)
                if len(obs_var) > 2:
                    continue
                if var_from_sanger[k].has_key(ind):  
                    var_from_sanger[k][ind].append(obs_var)
                else:
                    var_from_sanger[k][ind] = [obs_var]
                var_from_sanger[k]['REF'] = refseq.replace('-','')[pos-1]

    invar_pos = defaultdict(list)
    for ind,target,refstart,reflen,refseq,sangseq,m in sanger_data:
        pos = 0
        for rbase,sbase in zip(refseq,sangseq):
            if rbase != '-':
                pos += 1
            k = (target,str(refstart+pos))
            if k in var_from_sanger.keys():
                if not ind in var_from_sanger[k].keys():
                    obs_var = Seq.undegenerate_iupac(sbase)
                    if len(obs_var) <= 2:
                        var_from_sanger[k][ind] = [obs_var]
            elif rbase == sbase and rbase.upper() != 'N':
                invar_pos[k].append(ind)

    invar_pos_count = {}
    for k in invar_pos:
        invar_pos_count[k] = Util.countdict(invar_pos[k])

    return var_from_sanger, invar_pos, invar_pos_count

def get_metrics(vcf,add_QUAL=True):
    metrics = {}
    if add_QUAL:
        metrics['QUAL'] = ('Float','site quality')
    fh = open(vcf)
    for l in fh:
        if l.startswith('#CHROM'):
            return metrics
        elif l.startswith('##INFO'):
            m = re.search('ID=(.+?),Number=1,Type=(Integer|Float),Description="(.+?)"',l)
            if m:
                metrics[m.groups()[0]] = m.groups()[1:]

def logif(x):
    '''log10 '''
    return x<0 and -numpy.log10(-x) or (x and numpy.log10(x) or x)

def passfn(x):
    ''' '''
    return x

def gen_model_zero(test_metrics,model_terms):
    mod_zero = []
    for mt in test_metrics:
        mod_zero.append(model_terms[mt] == 'above' and -numpy.inf or numpy.inf)
    return mod_zero

def choose_bin(val,bins):
    for b in bins:
        if val <= b:
            return b
    return b


import time
def roc_calc(vcf, var_sites, invar_sites, model, values, verbose = True,eval_sites=None):
    if eval_sites is None:
        eval_sites = set(invar_sites+var_sites).intersection(set(vcf.keys()))
    X = []
    Y = []
    tot_var = len(var_sites)
    start_sec = time.time()
    for i,val in enumerate(values):
        if i % 10 == 0:
            elapsed = time.time() - start_sec
            pct_done = float(i) / len(values)
            if pct_done:
                est_remain = elapsed/pct_done - elapsed
            else:
                est_remain = numpy.inf
            if verbose: print >> sys.stderr, '\r%s / %s %0.1f %% %0.1f elapsed %0.1f est. remain' % (i,len(values),pct_done*100,elapsed,est_remain),
        pos_sites = [k for k in eval_sites if model(vcf[k],val)]
        try:
            X.append(float(len(set(pos_sites).intersection(set(invar_sites))))/len(pos_sites))
            Y.append(float(len(set(pos_sites).intersection(set(var_sites))))/tot_var)
        except:
            pass
    return X,Y

def hypotenuse(p1,p2):
    if not (len(p1) == 2 and len(p2) == 2):
        print >> sys.stderr, '%s and %s must be (x,y) points' % (p1,p2)
        raise TypeError, 'invalid arguments'
    import math
    x = p1[0] - p2[0]
    y = p1[1] - p2[1]
    return math.sqrt(x**2 + y**2)


def filter_roc(X,Y,test_values,min_dist=0.05):
    Yfilt = []
    Xfilt = []
    test_values_filt = []
    for x,y,cut in zip(X,Y,test_values):
        if len(Xfilt) == 0 or hypotenuse((Xfilt[-1],Yfilt[-1]),(x,y)) >= min_dist:
            Xfilt.append(x)
            Yfilt.append(y)
            test_values_filt.append(cut)

    return Xfilt,Yfilt,test_values_filt

def model_gen(metric,test):
    if test == 'above':
        return lambda sd,x: sd.has_key(metric) and float(sd[metric]) > x
    elif test == 'below':
        return lambda sd,x: sd.has_key(metric) and float(sd[metric]) < x
    elif test == 'absabove':
        return lambda sd,x: sd.has_key(metric) and abs(float(sd[metric])) > x
    elif test == 'absbelow':
        return lambda sd,x: sd.has_key(metric) and abs(float(sd[metric])) > x
    else:
        errstr = 'test must be in {above,below,absabove,absbelow}; given %s' % test
        raise ValueError, errstr

def calc_confusion_matrix_from_sanger(vcf_data, sanger_data, snpwin=5, positive_filter=lambda sd,ind:False, negative_filter=lambda sd,ind:False,skip_vcf_absences=True,use_indivs=None):
    '''
    given vcf per load_vcf and sanger per align_sanger, as well as two filters:

    positive_filter is applied to sites that differ from reference, negative_filter to sites that agree.  
    In either case the filters (functions) take a single value dict from vcf ("site dict" or "sd" per cutoff_fn of load_vcf)
    if return value is True, filters DO NOT permit the call (no call is made).  If return value is False, call is not filtered

    skips any variant in window with more than one variant (in sanger) in window size snpwin.
    '''
    xd = { '0':'REF','1':'ALT' }

    cm = {}.fromkeys(['TP','TN','FP','FN','P','N','Pfilt','Nfilt'],0)

    for ind,target,refstart,reflen,refseq,sangSeq,m in sanger_data:
        if use_indivs is not None and not ind in use_indivs:
            continue
        #print ind,target,refstart,reflen
        ngaps = 0
        for i in range(reflen):
            if refseq[i] == '-':
                ngaps += 1
                continue
            if refseq[i] == 'N' or sangSeq[i] == 'N':
                continue
            if len(m[i-(snpwin/2):i+(snpwin/2)+1].replace('.','')) > 1:
                #print '-',
                continue
            if not skip_vcf_absences:
                if m[i] == '.':
                    cm['N'] += 1
                else:
                    cm['P'] += 1
            if not vcf_data.has_key((target,str(i+refstart+1-ngaps))):
                #print (target,str(i+refstart+1-ngaps))
                continue
            if not vcf_data[(target,str(i+refstart+1-ngaps))]['indiv_gt'].has_key(ind):
                continue
            try:
            #if 1:
                sd = vcf_data[(target,str(i+refstart+1-ngaps))]
                s = [sd[xd[gt]] for gt in sd['indiv_gt'][ind]['GT'].split('/')]
                qs = sd['indiv_gt'][ind]['GQ']
                qn = int(float(qs))

                if '1' in sd['indiv_gt'][ind]['GT']:
                    variant = True
                else:
                    variant = False
                    
                if skip_vcf_absences:
                    if m[i] == '.':
                        cm['N'] += 1
                    else:
                        cm['P'] += 1

                if variant:
                    if positive_filter(sd,ind): #called variant; filter as positive
                        cm['Pfilt'] += 1
                    else:
                        if set(s)==set(Seq.undegenerate_iupac(sangSeq[i])):
                            #print sangSeq[i],
                            cm['TP'] += 1
                        else:
                            #print set(s),set(Seq.undegenerate_iupac(sangSeq[i])),qn,sd.get('QD',None),sd['indiv_gt'][ind]['DP'],'FP'
                            cm['FP'] += 1
                    
                else:
                    if negative_filter(sd,ind): #called invariant; filter as negative
                        cm['Nfilt'] += 1
                    else:
                        if set(s)==set(Seq.undegenerate_iupac(sangSeq[i])):
                            #print '.',
                            cm['TN'] += 1
                        else:
                            #print set(s),set(Seq.undegenerate_iupac(sangSeq[i])),qn,sd.get('QD',None),sd['indiv_gt'][ind]['DP'],'FN'
                            cm['FN'] += 1
            except:
            #else:
                pass
        #print '\n'

    return cm
    


def calc_confusion_matrix_from_sanger2(vcf_data, sanger_data, snpwin=5, \
                                       site_hi_bounds={},site_lo_bounds={},ind_hi_bounds={},ind_lo_bounds={},site_neg_param=[],ind_neg_param=['GQ','DP'],\
                                       skip_vcf_absences=True,use_indivs=None):
    '''
    given vcf per load_vcf and sanger per align_sanger, as well as two filters:

    positive_filter is applied to sites that differ from reference, negative_filter to sites that agree.  
    In either case the filters (functions) take a single value dict from vcf ("site dict" or "sd" per cutoff_fn of load_vcf)
    if return value is True, filters DO NOT permit the call (no call is made).  If return value is False, call is not filtered

    skips any variant in window with more than one variant (in sanger) in window size snpwin.
    '''
    xd = { '0':'REF','1':'ALT' }

    cm = {}.fromkeys(['TP','TN','FP','FN','P','N','Pfilt','Nfilt'],0)

    for ind,target,refstart,reflen,refseq,sangSeq,m in sanger_data:
        if use_indivs is not None and not ind in use_indivs:
            continue
        #print ind,target,refstart,reflen
        ngaps = 0
        for i in range(reflen):
            if refseq[i] == '-':
                ngaps += 1
                continue
            if refseq[i] == 'N' or sangSeq[i] == 'N':
                continue
            if len(m[i-(snpwin/2):i+(snpwin/2)+1].replace('.','')) > 1:
                #print '-',
                continue
            if not skip_vcf_absences:
                if m[i] == '.':
                    cm['N'] += 1
                else:
                    cm['P'] += 1
            if not vcf_data.has_key((target,str(i+refstart+1-ngaps))):
                #print (target,str(i+refstart+1-ngaps))
                continue
            if not vcf_data[(target,str(i+refstart+1-ngaps))]['indiv_gt'].has_key(ind):
                continue
            try:
            #if 1:
                sd = vcf_data[(target,str(i+refstart+1-ngaps))]
                s = [sd[xd[gt]] for gt in sd['indiv_gt'][ind]['GT'].split('/')]
                qs = sd['indiv_gt'][ind]['GQ']
                qn = int(float(qs))

                if '1' in sd['indiv_gt'][ind]['GT']:
                    variant = True
                else:
                    variant = False

                if filt_gt_call(sd,ind,site_hi_bounds,site_lo_bounds,ind_hi_bounds,ind_lo_bounds,site_neg_param,ind_neg_param):
                    filtered = True
                else:
                    filtered = False
                
                if skip_vcf_absences:
                    if m[i] == '.':
                        cm['N'] += 1
                    else:
                        cm['P'] += 1

                if variant:
                    if filtered:
                        cm['Pfilt'] += 1
                    else:
                        if set(s)==set(Seq.undegenerate_iupac(sangSeq[i])):
                            #print sangSeq[i],
                            cm['TP'] += 1
                        else:
                            #print set(s),set(Seq.undegenerate_iupac(sangSeq[i])),qn,sd.get('QD',None),sd['indiv_gt'][ind]['DP'],'FP'
                            cm['FP'] += 1
                    
                else:
                    if filtered:
                        cm['Nfilt'] += 1
                    else:
                        if set(s)==set(Seq.undegenerate_iupac(sangSeq[i])):
                            #print '.',
                            cm['TN'] += 1
                        else:
                            #print set(s),set(Seq.undegenerate_iupac(sangSeq[i])),qn,sd.get('QD',None),sd['indiv_gt'][ind]['DP'],'FN'
                            cm['FN'] += 1
            except:
            #else:
                pass
        #print '\n'

    return cm

def vcf_filter_fns(gq=0,qd=0,q=0,dp=0,fh=0):
    
    def pos_filt(sd,ind):
        return float(sd['indiv_gt'].get(ind,{}).get('GQ',0)) < gq \
               or int(sd['indiv_gt'].get(ind,{}).get('DP',0)) < dp \
               or float(sd.get('QD',0)) < qd \
               or float(sd.get('QUAL',0)) < q \
               or fract_het(sd['indiv_gt']) > fh
    def neg_filt(sd,ind):
        return float(sd['indiv_gt'].get(ind,{}).get('GQ',0)) < gq \
               or int(sd['indiv_gt'].get(ind,{}).get('DP',0)) < dp

    return pos_filt,neg_filt

def rates_from_confusion_matrix(cm):
    return float(cm['FP'])/(cm['N']),float(cm['TP'])/(cm['P'])
    #return float(cm['FP'])/(cm['TN']+cm['FP']),float(cm['TP'])/(cm['TP']+cm['FN'])

def filt_gt_call(sd,ind,site_hi_bounds={},site_lo_bounds={},ind_hi_bounds={},ind_lo_bounds={},site_neg_param=[],ind_neg_param=['GQ','DP']):
    
    if '1' in sd['indiv_gt'][ind]['GT']:
        site_param = list(set(site_hi_bounds.keys() + site_lo_bounds.keys()))
        ind_param = list(set(ind_hi_bounds.keys() + ind_lo_bounds.keys()))
        return any([not( site_lo_bounds.get(param,0) <= float(sd[param]) < site_hi_bounds.get(param,numpy.inf) ) for param in site_param]) \
               or any([not( ind_lo_bounds.get(param,0) <= float(sd['indiv_gt'][ind][param]) < ind_hi_bounds.get(param,numpy.inf) ) for param in ind_param])
    else:
        return any([not( site_lo_bounds.get(param,0) <= float(sd[param]) < site_hi_bounds.get(param,numpy.inf) ) for param in site_neg_param]) \
               or any([not( ind_lo_bounds.get(param,0) <= float(sd['indiv_gt'][ind][param]) < ind_hi_bounds.get(param,numpy.inf) ) for param in ind_neg_param])

    
def roc_helper(vcf,fas,ref,wins,GQs,QDs,QUALs,DPs,FHs,prune_points=True,point_spacing = 0.01,skip_vcf_absences=True,use_indivs=None,exp_coord_cut_fn=None,ind_lookup_fn=lambda x:x):
    '''
    FHs == fract het(0-1)
    '''
    from video_analysis.vidtools import hypotenuse
    rocs = {}
    conditions = [(w,gq,qd,q,dp,fh) for w in wins for gq in GQs for qd in QDs for q in QUALs for dp in DPs for fh in FHs]
    invarstr  = ''
    if len(wins) == 1:
        invarstr += '-w%s' % wins[0]
    if len(GQs) == 1:
        invarstr += '-gq%s' % GQs[0]
    if len(QDs) == 1:
        invarstr += '-qd%s' % QDs[0]
    if len(QUALs) == 1:
        invarstr += '-q%s' % QUALs[0]
    if len(DPs) == 1:
        invarstr += '-dp%s' % DPs[0]
    if len(FHs) == 1:
        invarstr += '-fh%s' % FHs[0]
    
    i = 0
    tickon = (len(conditions) * len(fas)) / 100
    if tickon < 1:
        tickon = 1
    print >> sys.stderr, '%s conditions by %s test sets = %s tests; dots on %s' % (len(conditions),len(fas),(len(conditions) * len(fas)), tickon)
    for fa in fas:
        sanger_data = align_sanger(fa,ref,exp_coord_cut_fn,ind_lookup_fn)
        p = []
        t = []
        for w,gq,qd,q,dp,fh in conditions:
            if i % tickon == 0: print >> sys.stderr, '.',
            i += 1
            def pos_filt(sd,ind):
                return float(sd['indiv_gt'].get(ind,{}).get('GQ',0)) < gq \
                       or int(sd['indiv_gt'].get(ind,{}).get('DP',0)) < dp \
                       or float(sd.get('QD',0)) < qd \
                       or float(sd.get('QUAL',0)) < q \
                       or fract_het(sd['indiv_gt']) > fh
            def neg_filt(sd,ind):
                return float(sd['indiv_gt'].get(ind,{}).get('GQ',0)) < gq \
                       or int(sd['indiv_gt'].get(ind,{}).get('DP',0)) < dp
            try:
                cm = calc_confusion_matrix_from_sanger( vcf,sanger_data,w,positive_filter=pos_filt,negative_filter=neg_filt,skip_vcf_absences=skip_vcf_absences,use_indivs=use_indivs)
                print >> sys.stderr, cm
                p.append( rates_from_confusion_matrix( cm ))
                varstr = ''
                if len(wins) != 1:
                    varstr += '-w%s' % w
                if len(GQs) != 1:
                    varstr += '-gq%s' % gq
                if len(QDs) != 1:
                    varstr += '-qd%s' % qd
                if len(QUALs) != 1:
                    varstr += '-q%s' % q
                if len(DPs) != 1:
                    varstr += '-dp%s' % dp
                if len(FHs) != 1:
                    varstr += '-fh%s' % fh
                
                
                t.append((p[-1][0],p[-1][1],"%s,tp%s,fp%s" % (varstr,cm['TP'],cm['FP'])))

            except:
                pass
        #print >> sys.stderr, p,t
        psort = [(0,0)]
        tsort = []
        #print >> sys.stderr, prune_points
        for pel,tel in sorted(zip(p,t)):
            #print >> sys.stderr,pel,tel,
            if prune_points == False or (pel[-1] > psort[-1][-1] and hypotenuse(pel,psort[-1]) > point_spacing):
                #print >> sys.stderr, 'added'
                psort.append(pel)
                tsort.append(tel)
        if len(fas) == 1:
            rocs[(invarstr,)] = (psort+[(1,1)],tsort)
        else:
            rocs[(fa,invarstr)] = (psort+[(1,1)],tsort)
    #x,y = Util.dezip(psort+[(1,1)])
    #plot(x,y,':',drawstyle="steps-post",label=fa)
    #null = [text(*te,**{'fontsize':6}) for te in tsort]
    return rocs

def plot_rocs(rocs,sym = '-',include_text=True,variable_as_label = False,label_prefix='',**kwargs):
    import Util
    from pylab import plot,text
    
    if variable_as_label:
        rocs_to_plot = {}
        allparam = []
        for k in rocs.keys():
            pli = re.findall('-([\w]+?)[\d\.]+',k[0])
            allparam.extend(pli)
            
        vparam = list(set(allparam))
        for k in rocs.keys():
            pli = re.findall('-([\w]+?)[\d\.]+',k[0])
            vp = [p for p in vparam if p not in pli]
            rocs_to_plot[(', '.join(vp),)] = rocs[k]
    else:
        rocs_to_plot=rocs


    
    for lbl,(psort,tsort) in rocs_to_plot.items():
        x,y = Util.dezip(psort+[(1,1)])
        if not kwargs.has_key('drawstyle'):
            kwargs['drawstyle'] = 'steps-post'
        plot(x,y,sym,lw=3,label=label_prefix+(''.join(lbl)),**kwargs)
        if include_text:
            null = [text(*te,**{'fontsize':6}) for te in tsort]

def write_structure_genotypes(vcf_data, outfile, keys_to_write = None, indiv_to_write = None):
    '''given a parsed vcf data structure per load_vcf and an outfile name

    prints a structure-formatted genotype file.

    if keys_to_write is supplied, only vcf_data items corresponding to those keys will be written
    '''

    if keys_to_write is None:
        keys_to_write = vcf_data.keys()

    keys_to_write.sort(key = lambda x: (x[0],int(x[1])))

    if indiv_to_write is None:
        indiv_to_write = set()
        for k in keys_to_write:
            v = vcf_data[k]
            indiv_to_write = indiv_to_write.union(set(v['indiv_gt'].keys()))
        indiv_to_write = sorted(list(indiv_to_write))

    ofh = open(outfile,'w')
    ofh.write('\t'.join(['%s.%s' % (c,p) for c,p in keys_to_write]))
    ofh.write('\n')

    for ind in indiv_to_write:
        ofh.write(ind)
        for k in keys_to_write:
            try:
                ofh.write('\t'+(vcf_data[k]['indiv_gt'][ind]['GT'].replace('/',' ')))
            except KeyError:
                ofh.write('\t-9 -9')
        ofh.write('\n')

    ofh.close()

        
def write_spagedi_genotypes(vcf_data, outfile, keys_to_write = None, indiv_to_write = None):
    '''generates output intended for SPAGeDi

    currently treats all individuals as originating from a single population;
    this will need to be elaborated upon
    '''

    from short_read_analysis import preprocess_radtag_lane
    lookup = dict([(l['sampleid'],l['population']) for l in preprocess_radtag_lane.get_table_as_dict('DB_library_data') if l.has_key('population')])

    if keys_to_write is None:
        keys_to_write = vcf_data.keys()

    keys_to_write.sort(key = lambda x: (x[0],int(x[1])))

    if indiv_to_write is None:
        indiv_to_write = set()
        for k in keys_to_write:
            v = vcf_data[k]
            indiv_to_write = indiv_to_write.union(set(v['indiv_gt'].keys()))
        indiv_to_write = sorted(list(indiv_to_write))

    ofh = open(outfile,'w')
    #write header
    ofh.write('%s\t1\t0\t%s\t1\t2\n0\nInd\tPop\t%s\n' % \
              (len(indiv_to_write),len(keys_to_write), '\t'.join(['%s.%s' % (c,p) for c,p in keys_to_write])))

    #write genotypes
    for ind in indiv_to_write:
        ofh.write('%s\t%s' % (ind,lookup.get(ind,'pop1')))
        for k in keys_to_write:
            try:
                gt = '/'.join([str(int(i)+1) for i in vcf_data[k]['indiv_gt'][ind]['GT'].split('/')])
                ofh.write('\t'+gt)
            except KeyError:
                ofh.write('\t0/0')
        ofh.write('\n')

    ofh.write('END\n')
    ofh.close()
    

def write_tassel_genotypes(vcf_data, outfile, keys_to_write = None, indiv_to_write = None):
    '''given vcf_data per load_vcf above and an outfile name

    writes tassel input format'''

    xd = { '0':'REF','1':'ALT' }

    if keys_to_write is None:
        keys_to_write = vcf_data.keys()

    keys_to_write.sort(key = lambda x: (x[0],int(x[1])))

    if indiv_to_write is None:
        indiv_to_write = set()
        for k in keys_to_write:
            v = vcf_data[k]
            indiv_to_write = indiv_to_write.union(set(v['indiv_gt'].keys()))
        indiv_to_write = sorted(list(indiv_to_write))

    ofh = open(outfile,'w')

    ofh.write('%s\t%s:%s\n' % (len(indiv_to_write),len(keys_to_write),'2'))
    if len(set([c for c,p in keys_to_write])) == 1:
        ofh.write('%s\n' % '\t'.join([p for c,p in keys_to_write]))
    else:
        ofh.write('%s\n' % '\t'.join(['%s.%s' % (c,p) for c,p in keys_to_write]))

    for ind in indiv_to_write:
        ofh.write('%s' % ind)
        for k in keys_to_write:
            try:
                gt = ':'.join([vcf_data[k][xd[i]] for i in vcf_data[k]['indiv_gt'][ind]['GT'].split('/')])
            except:
                gt = '?:?'
            ofh.write('\t' + gt)
        ofh.write('\n')

    ofh.close()
    
def write_plink_genotypes(vcf_data, outfile, keys_to_write = None, indiv_to_write = None):

    if keys_to_write is None:
        keys_to_write = vcf_data.keys()

    keys_to_write.sort(key = lambda x: (x[0],int(x[1])))

    if indiv_to_write is None:
        indiv_to_write = set()
        for k in keys_to_write:
            v = vcf_data[k]
            indiv_to_write = indiv_to_write.union(set(v['indiv_gt'].keys()))
        indiv_to_write = sorted(list(indiv_to_write))

    ofh = open(outfile+'.ped','w')

    idx = 0
    for ind in indiv_to_write:
        idx += 1
        ofh.write('FAM%s\t%s\t0\t0\t1\t0' % (idx,ind))
        for k in keys_to_write:
            try:
                gt = ' '.join([str(int(i)+1) for i in vcf_data[k]['indiv_gt'][ind]['GT'].split('/')])
            except:
                gt = '0 0'
            ofh.write('\t' + gt)
        ofh.write('\n')

    ofh.close()

    mapout = os.path.splitext(outfile)[0] + '.map'
    ofh = open(mapout,'w')
    infout = open(outfile + '.info','w')
    
    chrom_translation = dict([(k,i+1) for i,k in enumerate(set([k[0] for k in keys_to_write]))])
    open(mapout+'.xlat','w').write('\n'.join(['%s\t%s' % (k,v) for k, v in chrom_translation.items()]))

    for k in keys_to_write:
        ofh.write('%s\t%s.%s\t%s\t%s\n' % (chrom_translation[k[0]], k[0], k[1], 0, k[1]))
        infout.write('%s\t%s\n' % (k[1],k[1]))

    ofh.close()
    infout.close()


def write_flapjack_genotypes(vcf_data, outfile, keys_to_write = None, indiv_to_write = None):

    xd = { '0':'REF','1':'ALT' }


    if keys_to_write is None:
        keys_to_write = vcf_data.keys()

    keys_to_write.sort(key = lambda x: (x[0],int(x[1])))

    if indiv_to_write is None:
        indiv_to_write = set()
        for k in keys_to_write:
            v = vcf_data[k]
            indiv_to_write = indiv_to_write.union(set(v['indiv_gt'].keys()))
        indiv_to_write = sorted(list(indiv_to_write))

    ofh = open(outfile,'w')

    if len(set([c for c,p in keys_to_write])) == 1:
        ofh.write('\t%s\n' % '\t'.join([p for c,p in keys_to_write]))
    else:
        #ofh.write('\t%s\n' % '\t'.join(['%s.%s' % (c,p) for c,p in keys_to_write]))
        raise NotImplementedError, 'currently only supports single-chromosome output'

    for ind in indiv_to_write:
        ofh.write('%s' % ind)
        for k in keys_to_write:
            try:
                gt = '/'.join([vcf_data[k][xd[i]] for i in vcf_data[k]['indiv_gt'][ind]['GT'].split('/')])
            except:
                gt = '?/?'
            ofh.write('\t' + gt)
        ofh.write('\n')

    ofh.close()
    
    mapout = os.path.splitext(outfile)[0] + '.fj.map'
    ofh = open(mapout,'w')
    
    for k in keys_to_write:
        ofh.write('%s\t%s\t%s\n' % (k[1], 1, k[1]))

    ofh.close()

def write_peas_genotypes(vcf_data, outfile, keys_to_write = None, indiv_to_write = None):
    '''given vcf_data per load_vcf above and an outfile name

    writes peas input format'''

    xd = { '0':'REF','1':'ALT' }

    if keys_to_write is None:
        keys_to_write = vcf_data.keys()

    keys_to_write.sort(key = lambda x: (x[0],int(x[1])))

    if indiv_to_write is None:
        indiv_to_write = set()
        for k in keys_to_write:
            v = vcf_data[k]
            indiv_to_write = indiv_to_write.union(set(v['indiv_gt'].keys()))
        indiv_to_write = sorted(list(indiv_to_write))

    ofh = open(outfile,'w')

    ofh.write('SNPID\t%s\n' % '\t'.join('%s.%s' % (k[0],k[1]) for k in  keys_to_write))
    ofh.write('Chrom\t%s\n' % '\t'.join([k[0] for k in keys_to_write]))
    ofh.write('Position\t%s\n' % '\t'.join([k[1] for k in keys_to_write]))
    ofh.write('AlleleState\t%s\n' % '\t'.join(['%s/%s' % (vcf_data[k]['REF'],vcf_data[k]['ALT']) for k in keys_to_write]))
    ofh.write('Strand\t%s\n' % '\t'.join(['+'] * len(keys_to_write)))
    
    for ind in indiv_to_write:
        ofh.write('%s' % ind)
        for k in keys_to_write:
            try:
                gt = ''.join([vcf_data[k][xd[i]] for i in vcf_data[k]['indiv_gt'][ind]['GT'].split('/')])
            except:
                gt = '??'
            ofh.write('\t' + gt)
        ofh.write('\n')

    ofh.close()

def write_bimbam_genotypes(vcf_data, outfile, keys_to_write = None, indiv_to_write = None):
    '''does not suppport phased or multiple chromosome input vcfs
    '''

    xd = { '0':'REF','1':'ALT' }


    if keys_to_write is None:
        keys_to_write = vcf_data.keys()

    keys_to_write.sort(key = lambda x: (x[0],int(x[1])))

    if indiv_to_write is None:
        indiv_to_write = set()
        for k in keys_to_write:
            v = vcf_data[k]
            indiv_to_write = indiv_to_write.union(set(v['indiv_gt'].keys()))
        indiv_to_write = sorted(list(indiv_to_write))

    ofh = open(outfile,'w')

    if len(set([c for c,p in keys_to_write])) == 1:
        ofh.write('%s\n%s\n' % (len(indiv_to_write),len(keys_to_write)))
        ofh.write('IND,%s\n' % ','.join([ind for ind in indiv_to_write]))
    else:
        #ofh.write('\t%s\n' % '\t'.join(['%s.%s' % (c,p) for c,p in keys_to_write]))
        raise NotImplementedError, 'currently only supports single-chromosome output'

    for k in keys_to_write:
        ofh.write('%s' % k[1])
        for ind in indiv_to_write:
            try:
                gt = ''.join([vcf_data[k][xd[i]] for i in vcf_data[k]['indiv_gt'][ind]['GT'].split('/')])
            except:
                gt = '??'
            ofh.write(',' + gt)
        ofh.write('\n')

    ofh.close()
    
    mapout = os.path.splitext(outfile)[0] + '.bb.map.txt'
    ofh = open(mapout,'w')
    
    for k in keys_to_write:
        ofh.write('%s,%s\n' % (k[1], k[1]))

    ofh.close()

def write_fastPHASE_genotypes(vcf_data,outbase, keys_to_write = None, indiv_to_write = None):
    '''
    given vcf per load_vcf, writes one output file PER CHROM, name as:
    <outbase>_<CHROM>.inp

    it is recommended that outbase include path information to a directory in which these will be written, i.e.:

    outbase="beachmouse_fastPHASE/beaches_l55_k4_QD8_GQ12"

    will produce files:

    beachmouse_fastPHASE/beaches_l55_k4_QD8_GQ12_agouti_bac.inp
    beachmouse_fastPHASE/beaches_l55_k4_QD8_GQ12_mc1r_bac.inp
    ...etc
    
    directories will be created as necessary.
    '''

    
    xd = { '0':'REF','1':'ALT' }


    if keys_to_write is None:
        keys_to_write = vcf_data.keys()

    keys_to_write.sort(key = lambda x: (x[0],int(x[1])))

    if indiv_to_write is None:
        indiv_to_write = set()
        for k in keys_to_write:
            v = vcf_data[k]
            indiv_to_write = indiv_to_write.union(set(v['indiv_gt'].keys()))
        indiv_to_write = sorted(list(indiv_to_write))

    outroot = os.path.dirname(outbase)

    try:
        os.makedirs(outroot)
    except:
        pass


    this_chrom = None
    this_len = 0
    poslist = []
    for k in keys_to_write:
        if k[0] != this_chrom:
            if this_chrom is not None:
                outfile = '%s_%s.inp' % (outbase,this_chrom)
                posfile = '%s_%s.pos' % (outbase,this_chrom)
                open(posfile,'w').write(('\t'.join(poslist)) + '\n')
                ofh = open(outfile,'w')
                ofh.write('%s\n%s\n' % (len(indiv_to_write), this_len))
                for ind in indiv_to_write:
                    ofh.write('# %s\n' % ind)
                    ofh.write('%s\n%s\n' % tuple([''.join(li) for li in Util.dezip(chrom_data[ind])]))
                ofh.close()
            chrom_data = defaultdict(list)
            this_chrom = k[0]
            this_len = 0
            poslist = []
        for ind in indiv_to_write:
            try:
                chrom_data[ind].append(vcf_data[k]['indiv_gt'][ind]['GT'].split('/'))
            except:
                chrom_data[ind].append(['?','?'])
        this_len += 1
        poslist.append(k[1])

    #flush:
    outfile = '%s_%s.inp' % (outbase,this_chrom)
    posfile = '%s_%s.pos' % (outbase,this_chrom)
    open(posfile,'w').write(('\t'.join(poslist)) + '\n')
    ofh = open(outfile,'w')
    ofh.write('%s\n%s\n' % (len(indiv_to_write), this_len))
    for ind in indiv_to_write:
        ofh.write('# %s\n' % ind)
        ofh.write('%s\n%s\n' % tuple([''.join(li) for li in Util.dezip(chrom_data[ind])]))
    ofh.close()

from short_read_analysis import extract_genotypes_from_mclgr
def write_smartpca_genotypes(vcf_data,outbase):
    pm,gt = extract_genotypes_from_mclgr.genotypes_from_vcf_obj(vcf_data,0)
    extract_genotypes_from_mclgr.output_genotype_file(pm,gt,outbase+'.snp',outbase+'.ancestrymapgeno')
    open(outbase+'.ind','w').write('\n'.join(['%s\tU\t1' % ind for ind in gt.keys()]))
    open(outbase+'.par','w').write('genotypename:\t%s.ancestrymapgeno\n' \
                                   'snpname:\t%s.snp\n' \
                                   'indivname:\t%s.ind\n' \
                                   'evecoutname:\t%s.evec\n' \
                                   'evaloutname:\t%s.eval\n' \
                                   'snpweightoutname:\t%s.snpweightout\n' % tuple([outbase]*6))
