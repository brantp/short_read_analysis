#!/usr/bin/env python
'''given an mcl graph file and supporting information, uses "dirtiness" and fraction of individuals present to extract "well behaved" clusters and output genotypes for all individuals at those loci in the specified format
'''
import os,sys,re, Aln, Util
from collections import defaultdict
from Bio import pairwise2
from subprocess import Popen, PIPE
from cStringIO import StringIO

def load_cluster_data(gr,tab,mID=None):
    '''given an mcl graph output file, an mcl tab file of sequence labels, and an mID file of individuals per label,

    returns (clusters,labels,mID_by_label) dicts
    '''

    #extract clusters from mcl output
    print >>sys.stderr, 'load graph...',
    fc = open(gr).read()
    print >> sys.stderr, 'done\nparse graph...',
    body = re.search('begin(.+?)\)',fc.replace('\n',' ')).groups()[0]
    clusters = dict([(s.strip().split()[0],s.strip().split()[1:]) for s in body.strip().split('$') if len(s.strip().split()) > 1])
    print >> sys.stderr, 'done\nload labels...',

    #cluster labels
    labels = dict([l.strip().split() for l in open(tab).readlines()])

    print >> sys.stderr, 'done\nload mIDs...',
    #individuals by labels
    if mID is None:
        mID_by_label = None
    else:
        mID_by_label = dict([(l.split()[0],l.split()[1:]) for l in open(mID).readlines()])
    print >> sys.stderr, 'done'

    return clusters,labels,mID_by_label

def calculate_graph_properties(clusters,labels,mID_by_label,min_cov=2):
    '''given returns from load_cluster_data,

    prunes clusters of sequences with coverage (reads/individuals) < min_cov

    return (clusters,reads_per_cluster,indiv_per_cluster,cluster_dirt)

    calculates/returns:
    pruned clusters
    total reads per cluster
    total individuals per cluster
    "dirtiness" score:
        total reads for all indivs beyond top two sequences for EACH indiv / total reads per cluster
        e.g. if:
            l1.2: ["m1","m2"]
            l2.34: ["m2","m3"]
            l3.60: ["m2","m1"]
            dirt = 1 / 96
        "m1" considered het, "m3" homozygous, "m2" in violation by l1

        known issue: arguably overcounts dirtiness;
        takes as problematic the entire set of reads l1, not just the portion contributed by "m2"
        
    '''

    reads_per_cluster = {}
    indiv_per_cluster = {}
    cluster_dirt = {}
    numclust = len(clusters)
    tickon = numclust/100 + 1
    print >> sys.stderr, 'calculate read count properties for %s clusters' % numclust

    for i,cl in enumerate(clusters.keys()):
        if i%tickon==0: print >> sys.stderr, '\t%10d / %10d' % (i,numclust)

        #print >> sys.stderr, 'prune sequences below %s per-individual coverage' % min_cov
        keeps = []
        for l in clusters[cl]:
            lbl = labels[l]
            cnt = float(lbl.split('.')[1])
            if cnt/len(mID_by_label[lbl]) >= min_cov:
                keeps.append(l)
        if len(keeps) < 2:
            del(clusters[cl])
            continue
        else:
            clusters[cl] = keeps
        
        bad_labels = set([])
        all_mID = set(reduce(lambda x,y:x+y, [[m for m in mID_by_label[labels[l]] if not '-' in m] for l in clusters[cl]]))
        cnt_by_mID = {}
        for m in all_mID:
            cnt_by_mID[m] = []
            for l in clusters[cl]:
                lbl = labels[l]
                if m in mID_by_label[lbl]:
                    cnt = int(lbl.split('.')[1])
                    cnt_by_mID[m].append((cnt,lbl))

        #print >> sys.stderr, [labels[l] for l in clusters[cl]]
        #print >>sys.stderr, cnt_by_mID

        for m,cntli in cnt_by_mID.items():
            if len(cntli) > 2:
                for cnt,lbl in sorted(cntli,reverse=True)[2:]:
                    #print >> sys.stderr, '%s bad in %s' % (lbl,m)
                    bad_labels.add(lbl)
        #print >> sys.stderr, bad_labels, [int(lbl.split('.')[1]) for lbl in list(bad_labels)]

        dirt_num = sum([int(lbl.split('.')[1]) for lbl in list(bad_labels)])
        
        indiv_per_cluster[cl] = len(all_mID)
        reads_per_cluster[cl] = 0
        
        for l in clusters[cl]:
            lbl = labels[l]
            cnt = int(lbl.split('.')[1])
            reads_per_cluster[cl] += cnt

        cluster_dirt[cl] = dirt_num/float(reads_per_cluster[cl])

    return clusters,reads_per_cluster,indiv_per_cluster,cluster_dirt

def extract_alleles_from_cluster_labels(labels,snp_only=True,snp_flank=5):

    label_seqs = [(int(lbl.split('.')[1]),lbl,lbl.split('.')[2]) for lbl in labels]
    label_seqs = [l[1:] for l in sorted(label_seqs,reverse=True)]
    #print >> sys.stderr,'\n'.join([l[1] for l in label_seqs])

    mh = Popen(['muscle'],stdin=PIPE,stderr=PIPE,stdout=PIPE)
    mh.stdin.write('\n'.join(['>%s\n%s' % (l,s) for l,s in label_seqs]))
    alnstr = mh.communicate()[0]
    imfh = StringIO( alnstr )
    a = Aln.AlnFasta(imfh)
    #print alnstr
    #print a
    return a.extract_polymorphism_table(snp_only=snp_only,snp_flank=snp_flank)

def load_lines_from_uniqued(source_uniques,rv_sort = True, sort_key = lambda x: (len(x[0]),int(x[1])), keep_source_id = False):
    '''
    if keep_source_id is True
        returns list of 2-tuples uniqued_id (eg 100617_lane6_PE for "data/100617/100617_lane6_PE.uniqued")
        tuples are (parsed_lines,uniqued_id)

    else list of lines.
    '''
    uniquedlines = []
    for f in source_uniques:
        lines = []
        print >> sys.stderr, 'load %s ...' % f,
        lines = tuple([l.strip().split() for l in open(f).readlines()])
        print >> sys.stderr, '%s lines' % len(lines)
        if keep_source_id:
            uniqued_id = os.path.basename(os.path.splitext(f)[0])
            uniquedlines.extend( zip( lines,[uniqued_id]*len(lines) ) )
        else:
            uniquedlines.extend(lines)

    print >> sys.stderr, 'sort',
    if keep_source_id:
        uniquedlines.sort(reverse = rv_sort,key = lambda x: sort_key(x[0]))
    else:
        uniquedlines.sort(reverse = rv_sort,key = sort_key)
    print >> sys.stderr, 'done'
    return uniquedlines

def get_uniqued_lines_by_cluster(clid,clusters,labels,uniquedlines,remove_matched_lines=False):
    hits = []

    clseqs = [labels[lbl].split('.',3) for lbl in clusters[clid]]
    if remove_matched_lines:
        passed = []
        while uniquedlines:
            l = uniquedlines.pop(0)
            if any([l[0].startswith(seq[2]) for seq in clseqs]):
                hits.append(l)
            else:
                passed.append(l)
        
        return clseqs,hits,passed
    else:
        for l in uniquedlines:
            if any([l[0].startswith(seq[2]) for seq in clseqs]):
                hits.append(l)

        return clseqs,hits
            

def tabulate_genotypes(clusters, labels, mID_by_label, reads_per_cluster, indiv_per_cluster, cluster_dirt, indiv_cut = None, dirt_cut = 0.02):
    '''given a cutoff number of individuals which must be present to proceed with a marker,
    and the acceptable "dirtiness" (see calculate_graph_properties)

    returns polymorphism(=snp) (pm) and genotype (gt) dicts appropriate for output formatting as necessary
    pm keys are cluster ids, gt keys are individuals (inner dict: cluster ids)
    '''

    if indiv_cut is None:
        max_indiv = max(indiv_per_cluster.values())
        indiv_cut = int(0.9*max_indiv)

    pass_clust = [c for c in clusters.keys() if cluster_dirt.get(c,1) < dirt_cut and indiv_per_cluster.get(c,0) >= indiv_cut]

    print >> sys.stderr, 'clusters passing: %s with at least %s individuals and dirt cutoff of %s' % (len(pass_clust), indiv_cut, dirt_cut)

    ind = set([])
    for cl in pass_clust:
        for lbl in clusters[cl]:
            for m in mID_by_label[labels[lbl]]:
                ind.add(m)

    gt = {}
    for m in ind:
        gt[m] = defaultdict(dict)

    pm = {}

    tickon = len(pass_clust)/20

    for itval,cl in enumerate(pass_clust):
        if itval % tickon == 0: print >> sys.stderr, '%5d / %5d' % (itval,len(pass_clust))
        gt_by_mID = defaultdict(list)
        for lbk in clusters[cl]:
            lbl = labels[lbk]
            if len(mID_by_label[lbl]) > 1:
                for m in mID_by_label[lbl]:
                    gt_by_mID[m].append(lbl)
        #sorts within each individual's list of recovered haplotypes by the number of times the haplotype was observed OVERALL
        [v.sort(key = lambda x: int(x.split('.')[1]), reverse=True) for v in gt_by_mID.values()]

        haplotype_poly_table = extract_alleles_from_cluster_labels([labels[lbk] for lbk in clusters[cl]])

        for i,pos in enumerate(haplotype_poly_table[0]):
            alleles = defaultdict(int)
            site_gt_by_mID = {}
            for m,v in gt_by_mID.items():
                gt_raw = [haplotype_poly_table[1][lbl][i] for lbl in v]
                these_gt = sorted(list(set(gt_raw)),key=lambda x: gt_raw.index(x))
                
                if len(these_gt) == 1:
                    alleles[these_gt[0]] += 2
                    site_gt_by_mID[m] = these_gt*2
                else:
                    for al in these_gt[:2]:
                        alleles[al] += 1
                    site_gt_by_mID[m] = these_gt[:2]
                    polarized_loci,polarized_geno = extract_genotypes_from_mclgr.genotypes_by_parent(pm,gt,parents,hybrids=hybrids,remove_targets=reduce(lambda x,y: x+y,parents.values()) + hybrids)
            if len(alleles) == 1:
                continue

            ranked_alleles = [k for k,v in sorted(alleles.items(),key = lambda x:x[1],reverse=True)]

            site_name = '%s.%s' % (cl,haplotype_poly_table[0][i])
            pm[site_name] = ranked_alleles
            for m,gtl in sorted(site_gt_by_mID.items()):
                this_gt = [ranked_alleles.index(g) for g in gtl]
                gt[m][site_name] = this_gt

            

            ##########
            # display only
            '''
            print >> sys.stderr, '\n','%s.%s' % (cl,haplotype_poly_table[0][i]),sorted([(v,k) for k,v in alleles.items()],reverse=True),'\n',ranked_alleles

            count_gt = defaultdict(int)
            count_gt = defaultdict(int)
            for m,gtl in sorted(site_gt_by_mID.items()):
                this_gt = [ranked_alleles.index(g) for g in gtl]
                print >> sys.stderr,m,this_gt
                count_gt[m[:2]+str(sorted(this_gt))]+=1
            print >> sys.stderr,sorted(count_gt.items())
            '''
            #
            #########

    return pm,gt        


def output_genotype_file(pm,gt,snpfn,gtfn):
    '''

    generates ancestrymapgeno and snp files for smartpca
    
    use something like
    for m in gt.keys(): 
        for s in gt[m].keys(): 
            print >> fh, '%s\t%s\t%s' % (s,m,int(sum([1 for a in gt[m][s] if a])))

    for id in pm.keys(): 
        print >> fh, '%s\t1\t0.0\t0\t%s' % (id, '\t'.join([str(c) for c in pm[id]]))

    '''
    fh = open(gtfn,'w')
    for m in gt.keys(): 
        for s in gt[m].keys(): 
            print >> fh, '%s\t%s\t%s' % (s,m,int(sum([1 for a in gt[m][s] if a == 0])))

    fh = open(snpfn,'w')
    for id in pm.keys(): 
        print >> fh, '%s\t1\t0.0\t0\t%s' % (id, '\t'.join([str(c) for c in pm[id]]))


def genotypes_from_vcf_obj(vcf,min_indiv=40):
    pm = {}
    gt = defaultdict(dict)
    for k,v in vcf.items():
        if len(v['indiv_gt']) < min_indiv or v['mac'] == 0: continue
        sname = '%s.%s' % k
        pm[sname] = [v['REF'],v['ALT']]
        for ind,gtdict in v['indiv_gt'].items():
            gt[ind][sname] = [int(i) for i in gtdict['GT'].split('/')]

    return pm,gt


def genotypes_by_parent(pm,gt,parents,hybrids=None,remove_targets=None):
    '''Given pm and gt per tabulate_genotypes,
    and dictionary parents of the form
    {'A': ['BW1','BW2']...}
    + optional list hybrids
    returns loci and genotypes suitable for output_cross_radtag_genotypes
    + a list of target individuals to remove (e.g. parents, problematic indivs)'''

    print >> sys.stderr, '%s snps to evalute' % (len(pm))

    fixedsnps = {}
    for site in pm.keys():
        gts = defaultdict(list)
        for pk in parents.keys():
            for indiv in parents[pk]:
                try:
                    gts[pk].extend(gt[indiv].get(site,[]))
                except:
                    print indiv,site
                    raise
        gts['A'] = list(set(gts['A']))
        gts['B'] = list(set(gts['B']))
        if len(gts['A']) == 1 and len(gts['B']) == 1 and gts['A'][0] != gts['B'][0]:
            gtlookup = {}
            for pk in gts.keys():
                gtlookup[gts[pk][0]] = pk
            fixedsnps[site] = gtlookup
    print >> sys.stderr, '%s fixed snps between parents' % (len(fixedsnps))
    
    if hybrids:
        f1het = []
        for site in fixedsnps.keys():
            #print >>sys.stderr,[gt[polarized_loci,polarized_geno = extract_genotypes_from_mclgr.genotypes_by_parent(pm,gt,parents,hybrids=hybrids,remove_targets=reduce(lambda x,y: x+y,parents.values()) + hybrids)h][site] for h in hybrids if gt[h][site]]
            #print >>sys.stderr,[(len(gt[h][site]), set([fixedsnps[site][gt[h][site][i]] for i in range(2)])) for h in hybrids if gt[h][site]]
            if all([len(gt[h].get(site,[]))==2 and set(['A','B'])==set([fixedsnps[site][gt[h].get(site,[])[i]] for i in range(2)]) for h in hybrids if gt[h].get(site,[])]):
                f1het.append(site)
        loci = f1het
        print >> sys.stderr, '%s parent-fixed snps het in hybrids' % (len(f1het))        
    else:
        loci = fixedsnps.keys()


    genotypes = {}
    for indiv in gt.keys():
        genotypes[indiv] = {}
        for site in fixedsnps.keys():
            g = ''.join(sorted([fixedsnps[site].get(allele,'') for allele in gt[indiv].get(site,[])]))
            if len(g) == 2:
                genotypes[indiv][site] = g
    if remove_targets is not None:
        for t in remove_targets:
            del genotypes[t]

    return loci, genotypes

def get_sex_info():
    from short_read_analysis import preprocess_radtag_lane
    db = preprocess_radtag_lane.get_table_as_dict('DB_library_data')
    sex = dict([(l['sampleid'],l['sex']) for l in db if l.has_key('sex')])
    hyb_parent_sex = {}
    for l in db:
        if l.has_key('sire') and l.has_key('dam'):                      
            if l['sire'].startswith('F1') and l['dam'].startswith('BW'):
                hyb_parent_sex[l['sampleid']] = 'M'
            elif l['dam'].startswith('F1') and l['sire'].startswith('BW'):
                hyb_parent_sex[l['sampleid']] = 'F'

    return hyb_parent_sex, sex
            
gt_classes = {'AA':0,'AB':1,'BB':2}

def get_bc_obs_by_parent(site,polarized_geno,sex,hyb_parent_sex):
    obs_mat = {}
    for ps in ['M','F']:
        for s in ['M','F']:
            obs_mat[(ps,s)] = defaultdict(int)
    for ind in hyb_parent_sex.keys():
        if polarized_geno.get(ind,{}).has_key(site): obs_mat[(hyb_parent_sex[ind],sex[ind])][gt_classes[polarized_geno[ind][site]]] += 1
    return obs_mat

def get_bc_X_exp_by_parent(obs_mat):
    exp_mat = {}
    exp_mat[('F','F')] = {0: sum(obs_mat[('F','F')].values()) * 0.5, 1: sum(obs_mat[('F','F')].values()) * 0.5, 2: sum(obs_mat[('F','F')].values()) * 0.0}
    exp_mat[('F','M')] = {0: sum(obs_mat[('F','M')].values()) * 0.5, 1: sum(obs_mat[('F','M')].values()) * 0, 2: sum(obs_mat[('F','M')].values()) * 0.5}
    exp_mat[('M','F')] = {0: sum(obs_mat[('M','F')].values()) * 1, 1: sum(obs_mat[('M','F')].values()) * 0, 2: sum(obs_mat[('M','F')].values()) * 0}
    exp_mat[('M','M')] = {0: sum(obs_mat[('M','M')].values()) * 1, 1: sum(obs_mat[('M','M')].values()) * 0, 2: sum(obs_mat[('M','M')].values()) * 0}
    return exp_mat

def get_bc_aut_exp_by_parent(obs_mat):
    exp_mat = {}
    exp_mat[('F','F')] = {0: sum(obs_mat[('F','F')].values()) * 0.5, 1: sum(obs_mat[('F','F')].values()) * 0.5, 2: sum(obs_mat[('F','F')].values()) * 0.0}
    exp_mat[('F','M')] = {0: sum(obs_mat[('F','F')].values()) * 0.5, 1: sum(obs_mat[('F','F')].values()) * 0.5, 2: sum(obs_mat[('F','F')].values()) * 0.0}
    exp_mat[('M','F')] = {0: sum(obs_mat[('F','F')].values()) * 0.5, 1: sum(obs_mat[('F','F')].values()) * 0.5, 2: sum(obs_mat[('F','F')].values()) * 0.0}
    exp_mat[('M','M')] = {0: sum(obs_mat[('F','F')].values()) * 0.5, 1: sum(obs_mat[('F','F')].values()) * 0.5, 2: sum(obs_mat[('F','F')].values()) * 0.0}
    return exp_mat

def get_bc_y_exp_by_parent(obs_mat):
    exp_mat = {}
    exp_mat[('F','F')] = {0: sum(obs_mat[('F','F')].values()) * 0, 1: sum(obs_mat[('F','F')].values()) * 0, 2: sum(obs_mat[('F','F')].values()) * 0}
    exp_mat[('F','M')] = {0: sum(obs_mat[('F','F')].values()) * 1, 1: sum(obs_mat[('F','F')].values()) * 0, 2: sum(obs_mat[('F','F')].values()) * 0}
    exp_mat[('M','F')] = {0: sum(obs_mat[('F','F')].values()) * 0, 1: sum(obs_mat[('F','F')].values()) * 0, 2: sum(obs_mat[('F','F')].values()) * 0}
    exp_mat[('M','M')] = {0: sum(obs_mat[('F','F')].values()) * 0, 1: sum(obs_mat[('F','F')].values()) * 0, 2: sum(obs_mat[('F','F')].values()) * 1}
    return exp_mat


def add_pseudocount(cat_mat):
    for vout in cat_mat.values():
        for k in [0,1,2]:
            vout[k] += 1.0

def Xlinked_in_BC_chi2(site,polarized_geno,sex,hyb_parent_sex):
    obs_mat = get_bc_obs_by_parent(site,polarized_geno,sex,hyb_parent_sex)
    exp_mat = get_bc_X_exp_by_parent(obs_mat)
    add_pseudocount(obs_mat)
    add_pseudocount(exp_mat)
    o,e = [reduce(lambda x,y: x+y, l) for l in Util.dezip([([v for k,v in sorted(obs_mat[cat].items())],[v for k,v in sorted(exp_mat[cat].items())]) for cat in obs_mat.keys()])]
    lc2 = list_chi2(o,e)
    return sum(lc2)

def aut_in_BC_chi2(site,polarized_geno,sex,hyb_parent_sex):
    obs_mat = get_bc_obs_by_parent(site,polarized_geno,sex,hyb_parent_sex)
    exp_mat = get_bc_aut_exp_by_parent(obs_mat)
    add_pseudocount(obs_mat)
    add_pseudocount(exp_mat)
    o,e = [reduce(lambda x,y: x+y, l) for l in Util.dezip([([v for k,v in sorted(obs_mat[cat].items())],[v for k,v in sorted(exp_mat[cat].items())]) for cat in obs_mat.keys()])]
    lc2 = list_chi2(o,e)
    return sum(lc2)

def Ylinked_in_BC_chi2(site,polarized_geno,sex,hyb_parent_sex):
    obs_mat = get_bc_obs_by_parent(site,polarized_geno,sex,hyb_parent_sex)
    exp_mat = get_bc_y_exp_by_parent(obs_mat)
    add_pseudocount(obs_mat)
    add_pseudocount(exp_mat)
    o,e = [reduce(lambda x,y: x+y, l) for l in Util.dezip([([v for k,v in sorted(obs_mat[cat].items())],[v for k,v in sorted(exp_mat[cat].items())]) for cat in obs_mat.keys()])]
    lc2 = list_chi2(o,e)
    return sum(lc2)


def unit_chi2(o,e):

    if e > 0 :
        return (o-e)**2/e
    else:
        return 0

def list_chi2(obs,exp):
    '''given expected and observed counts lists

    returns chi2
    '''

    return [unit_chi2(o,e) for o,e in zip(obs,exp)]


def filter_Xlinked_loci(polarized_loci, polarized_geno, chi2_crit=20):
    hyb_parent_sex, sex = get_sex_info()
    xsites = []
    autsites = []
    for site in polarized_loci:
        if Xlinked_in_BC_chi2(site,polarized_geno,sex,hyb_parent_sex) < chi2_crit:
            xsites.append(site)
        else:
            autsites.append(site)

    return xsites,autsites
    
def output_cross_radtag_genotypes(loci,genotypes,filename,lg0='X'):
    '''Given list loci and dictionary genotype per genotypes_by_parent, writes file <filename>
    suitable for RQTL

    overloads 20101202:
    - if loci is a dict per maploci from load_cross_radtag_genotypes below, sort by map position in output
    - if filename is not string, use as filehandle (permits passing sys.stdout, for instance)
    '''

    def sortkey(x):
        if x == '':
            return 0
        else:
            return x

    if isinstance(loci,list):
        locnames = loci
        lgs = ['1' for i in range(len(loci))]
        mps = [str(i+1) for i in range(len(loci))]
    elif isinstance(loci,dict):
        for k,v in loci.items():
            if v[0] == 0:
                loci[k] = (lg0,v[1])
        locnames,lgs,mps = Util.dezip([(loc,str(lg),str(mp)) for loc,(lg,mp) in sorted(loci.items(),key=lambda x:[sortkey(v) for v in x[1]])])

    mID_lookup = dict([(m,str(i)) for i,m in enumerate(sorted(genotypes.keys()))])

    if isinstance(filename,str):
        fh = open(filename ,'w')
        #open(filename+'.mIDlookup','w').write('\n'.join(['%s\t%s' % (i,m) for m,i in sorted(mID_lookup.items())]))
    else:
        fh = filename
        #open(filename.name+'.mIDlookup','w').write('\n'.join(['%s\t%s' % (i,m) for m,i in sorted(mID_lookup.items())]))
        
    fh.write('ID,')
    fh.write(','.join(['%sr' % l for l in locnames]))
    fh.write('\n')
    fh.write(',')
    fh.write(','.join(lgs))
    fh.write('\n')
    fh.write(',')
    fh.write(','.join(mps))
    fh.write('\n')

    out_geno = {}
    
    for mID in genotypes.keys():
        #fh.write(mID_lookup[mID]+',')
        fh.write(mID+',')
        out_geno[mID] = dict([(mkr,genotypes[mID][mkr]) for mkr in locnames if genotypes[mID].has_key(mkr)])
        fh.write(','.join([genotypes[mID].get(mkr,'-') for mkr in locnames]))
        fh.write('\n')

    fh.close()

    return out_geno,mID_lookup


def load_cross_radtag_genotypes(gt_csv,mIDlookup=None,missing='-',id_header='ID'):
    '''given files of the form written by output_cross_radtag_genotypes,
    returns maploci and genotypes objects

    N.B. genotypes is identical to genotypes above,
    maploci is dict where keys are loci above,
    values are (linkage_group,genetic_position)'''

    fh = open(gt_csv,'U')

    if mIDlookup == 'skip':
        raise ValueError, '"skip" is deprecated; use False instead'
        #print >> sys.stderr, 'mIDlookup not used; skip'
        #mID_lookup = {}
        
    if mIDlookup is None:
        mIDlookup = gt_csv+'.mIDlookup'
        print >> sys.stderr, 'mIDlookup not specified, load from %s' % mIDlookup,
        if os.path.exists(mIDlookup):
            print >> sys.stderr,'successful'
        else:
            print >> sys.stderr,'FAILED'
            raise IOError, 'cannot find %s' % mIDlookup

    if mIDlookup:
        mID_lookup = dict([l.strip().split() for l in open(mIDlookup) if len(l.strip().split()) == 2 ])
        print >> sys.stderr, 'loaded mIDlookup data for %s individuals' % len(mIDlookup)
    else:
        mID_lookup = {}
    
    headrow = fh.readline().strip().split(',')
    chromrow = fh.readline().strip().split(',')
    posrow = fh.readline().strip().split(',')

    if id_header is None:
        lookup_idx = 0
    else:
        lookup_idx = headrow.index(id_header)

    print >> sys.stderr, 'using field %s of %s as ID' % (lookup_idx,len(headrow))

    #loci = fh.readline().strip().split(',')[1:]
    #maploci = dict(zip([l.rstrip('r') for l in loci],zip(map(int,fh.readline().strip().split(',')[1:]),map(float,fh.readline().strip().split(',')[1:]))))

    maploci = {}
    for header,chrom,pos in zip(headrow,chromrow,posrow):
        if len(chrom) > 0:
            if chrom == 'X':
                chrom = 0
            maploci[header.rstrip('r').lstrip('X')] = (int(chrom),float(pos))

    print >> sys.stderr, 'map data for %s sites loaded' % len(maploci)
    
    genotypes = defaultdict(dict)    

    for l in fh:
        fields = l.strip().split(',')
        for header,chrom,val in zip(headrow,chromrow,fields):
            if len(chrom) > 0 and val != missing:
                genotypes[mID_lookup.get(fields[lookup_idx],fields[lookup_idx])][header.rstrip('r').lstrip('X')] = val

    print >> sys.stderr, 'genotype data for %s individuals loaded' % len(genotypes)

    return maploci,genotypes
        
        

def calculate_coverage_by_site(loci,genotypes):
    '''calculate the number of individuals covered at a site as a dictionary'''

    coverage = {}
    for l in loci:
        coverage[l] = len(filter(None,[v.get(l,None) for v in genotypes.values()]))

    return coverage




def cross_allele_counts(loci,genotypes):
    '''measures allele frequencies and sorts based on a filter for
    minfreq and maxfreq of the A allele'''

    allele_counts = {}
    for l in loci:
        allele_counts[l] = [sum([v.get(l,'').count('A') for v in genotypes.values()]),sum([v.get(l,'').count('B') for v in genotypes.values()])]
        
    return allele_counts


def cross_alleles_passed(allele_counts,minfreq=None,maxfreq=None):

    alleles_passed = []
    for loc in allele_counts:
        if minfreq<=(allele_counts[loc][0]/float(allele_counts[loc][0]+allele_counts[loc][1]))<=maxfreq:
            alleles_passed.append(loc)

    return alleles_passed



def find_best_nonrepeated_snps(alleles_passed,coverage):
    '''if multiple snps at different positions in the same cluster,
    find_best_nonrepeated_snps'''

    keep = {}
    for site in alleles_passed:
        clust, pos = site.split('.')
        if clust in keep.keys():
            if keep[clust][0] < coverage[site]:
                keep[clust] = [coverage[site],site]
        else:
            keep[clust] = [coverage[site],site]

    return keep

    




if __name__ == "__main__":
    pass
