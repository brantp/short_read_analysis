'''functionality for enclosure experiments

prefix refers to string identifying indivs in sequencing samples (and therefore vcf files)

includes:
generate phenotype files and handle generation of genotypes (using variant_detection.write_wigs_output)
'''
try:
    from rtd.preprocess_radtag_lane import get_table_as_dict
except:
    from radtag_denovo.preprocess_radtag_lane import get_table_as_dict

import numpy,os,sys,re
from glob import glob

def indivs_in_vcf(vcf_obj):
    indivs = set([])
    for v in vcf_obj.values():
        these_indivs = set(v['indiv_gt'].keys())
        indivs = indivs.union(these_indivs)
    return indivs

def indiv_lists_by_enc(vcf_obj,db_name,query_d,prefix='RB',sort_key=None):
    '''
    see write_wigs_all_simple for example of query_d
    '''
    td = get_table_as_dict(db_name,suppress_fc_check=True)
    vcf_indivs = indivs_in_vcf(vcf_obj)
    if 'site' in query_d.keys():
        sites = query_d['site']
        indiv_lists = tuple([[d['id'] for d in td if d['site'] == this_site and prefix+d['id'] in vcf_indivs and all([d[k] in v for k,v in query_d.items() if k != 'site'])] for this_site in sites])
    else:
        indiv_lists = tuple([[d['id'] for d in td if prefix+d['id'] in vcf_indivs and all([d[k] in v for k,v in query_d.items()])]])
    if sort_key is not None:
        indiv_lists = [sorted(l,key=sort_key) for l in indiv_lists]
    return indiv_lists

def wigs_phenotypes_from_DB(db_name,indiv_lists,survive_field='recap1',survive_value='R'):
    '''indiv_lists is tuple of lists, one list per pop, of individual ids
    '''
    td = get_table_as_dict(db_name,suppress_fc_check=True)
    phenos = []
    for pop_n,indivs in enumerate(indiv_lists):
        phenos.append([])
        for ind_n,indiv in enumerate(indivs):
            phenos[-1].append([int(d[survive_field] == survive_value) for d in td if d['id'] == indiv][0])
    return phenos
            
def write_wigs_pheno(phenos,outfile):
    fh = open(outfile,'w')
    for pl in phenos:
        fh.write(' '.join(map(str,pl)))
        fh.write('\n')
    fh.close()

from variant_detection import write_wigs_genotypes

def write_wigs_all_simple(vcf_obj,db_name,query_d,outbase,by_locus=False,**kwargs):
    '''
    query_d dictionary represents columns and values from <db_name> spreadsheet to require for inclusion
    for example, query_d={'site':['L'],'enc':['SW']}
    '''
    indiv_lists = indiv_lists_by_enc(vcf_obj,db_name,query_d)
    phenos = wigs_phenotypes_from_DB(db_name,indiv_lists,**kwargs)

    open(outbase+'-wigs-indivs.tuple','w').write(indiv_lists.__repr__())
    write_wigs_pheno(phenos,outbase+'-wigs-pheno.txt')
    
    if by_locus:
        loci = list(set([c for c,s in vcf_obj.keys()]))
        outbase_by_locus = outbase+'-by_locus'
        print >> sys.stderr, 'write %s loci to %s' % (len(loci),outbase_by_locus)
        if not os.path.exists(outbase_by_locus): os.makedirs(outbase_by_locus)
        for i,loc in enumerate(loci):
            print >> sys.stderr, '\r %s / %s' % (i+1,len(loci)),
            this_vcf_obj = dict([((c,s),sd) for (c,s),sd in vcf_obj.items() if c == loc])
            this_outbase = os.path.join(outbase_by_locus,loc)
            write_wigs_genotypes(this_vcf_obj,indiv_lists,this_outbase,[0],[0],'RB')
    else:
        write_wigs_genotypes(vcf_obj,indiv_lists,outbase,[0],[0],'RB')

def cut_fn(sd): #NOTE QD 5 DIFFERS FROM SINERGIA/CRL
    summ_stats = ['FS','QUAL','BaseQRankSum','QD','SB','ReadPosRankSum']
    FS,QUAL,BaseQRankSum,QD,SB,ReadPosRankSum = map(float,[sd.get(ss,0) for ss in summ_stats])
    return FS < 72 and QUAL > 218 and BaseQRankSum < 5 and QD >= 5 and SB < 10 and ReadPosRankSum > -9 and sd['fh'] < 0.55

def parse_wigs(wigsout):
    premat = []
    for l in open(wigsout):
        premat.append(l.strip().split())
    mat = numpy.array(premat,dtype=float)
    return mat.transpose()

def sig_sites(wigsmat,lo=2.5,hi=97.5):
    sig = numpy.zeros(len(wigsmat),dtype=bool)
    for i,vec in enumerate(wigsmat):
        if numpy.mean(vec)>0:
            if numpy.percentile(vec,lo) > 0:
                sig[i] = True
        else:
            if numpy.percentile(vec,hi) < 0:
                sig[i] = True
    return sig

def num_sig_from_dir(indir,suf='s.txt',lo=2.5,hi=97.5,fdr=True):
    wigsouts = glob(os.path.join(indir,'*'+suf))
    print len(wigsouts)
    num_sig = sig_sites(parse_wigs(wigsouts[0]),lo,hi).astype(int)
    print num_sig.shape
    for i,wout in enumerate(wigsouts[1:]):
        try:
            num_sig += sig_sites(parse_wigs(wout),lo,hi).astype(int)
            print >> sys.stderr, '\r%s / %s' % (i+2,len(wigsouts)),
        except:
            print >> sys.stderr, '\r%s / %s FAILED' % (i+2,len(wigsouts))
    if fdr:
        return num_sig.astype(float) / float(len(wigsouts))
    else:
        return num_sig



def get_wigs_results(wigs_out,fdr_permutes_dir=None,site_coords=None,ofh=None):

    def write_summary(fh,coords,wS,sig,fdr=None):
        for i in range(len(wS)):
            l = '%s\t%s\t%s\t%s\t%s\n' % (coords[i][0],coords[i][1],wS[i],sig[i].astype(int),fdr is not None and fdr[i] or 0)
            fh.write(l)
        
            

    if site_coords is None:
        site_coords = re.sub('wigs_output.+?\.txt','wigs-sites.txt',wigs_out)
        if not os.path.exists(site_coords):
            print >> sys.stderr, 'attempted to use %s for coords; does not exist' % site_coords
            raise OSError
    coords = [eval(l) for l in open(site_coords).readlines()]

    print >> sys.stderr, 'load wigs output'
    wmat = parse_wigs(wigs_out)
    wS = wmat.mean(1)
    print >> sys.stderr, '%s wigs matrix loaded; %s sites mean: %s min: %s max: %s' % (wmat.shape,len(wS),wS.mean(),wS.min(),wS.max())
    print >> sys.stderr, 'calculate non-zero posteriors'
    sig = sig_sites(wmat)
    print >> sys.stderr, '%s sites non-zero (95%% credible)' % sig.sum()

    if fdr_permutes_dir is not None:
        print >> sys.stderr, 'load permutation results for FDR from %s' % fdr_permutes_dir 
        fdr = num_sig_from_dir(fdr_permutes_dir)
        nsig = (((fdr>0.01)+(1-sig))==0).sum()
        print >> sys.stderr, '1%% FDR sites: %s' % nsig
        if ofh is not None:
            write_summary(ofh,coords,wS,sig,fdr)
        return coords,wS,sig,fdr
    else:
        if ofh is not None:
            write_summary(ofh,coords,wS,sig)
        return coords,wS,sig

def compile_credible_sites_from_wigs_outputs(wigsouts,new_geno_base):
    '''wigsouts is list of wigs output-s files (e.g. from glob)

    new_geno base should be a filename base.
    new geno will add -credible-wigs-geno.txt;
    new sites will add -credible-wigs-sites.txt
    '''
    cred_coords = []
    geno_body = new_geno_base+'-credible-wigs-geno-body.txt'
    outfh = open(geno_body ,'w')
    header = ['','','','']
    for wo in wigsouts:
        geno = wo.replace('_output-s','-geno')
        coords,wS,sig = get_wigs_results(wo)
        idx_offset = len(cred_coords)
        cred_coords.extend(map(tuple,numpy.array(coords)[sig]))
        cred_sites = numpy.arange(len(sig))[sig]
        idx_lookup = dict(zip(cred_sites,range(len(cred_sites))))
        ifh = open(geno)
        for i in range(4):
            header[i] = ifh.next()
        for l in ifh:
            if int(l.split()[1]) in cred_sites:
                newl = l.split()
                newl[1] = str(idx_lookup[int(newl[1])] + idx_offset)
                outfh.write(' '.join(newl)+'\n')

    outfh.close()
    outfh = open(new_geno_base+'-credible-wigs-sites.txt','w')
    for coord in cred_coords:
        outfh.write(coord.__repr__()+'\n')
    outfh.close()

    geno_head = new_geno_base+'-credible-wigs-geno-head.txt'
    outfh = open(geno_head ,'w')
    outfh.write('%s %s\n' % (len(cred_coords),header[0].split()[1]))
    for i in range(1,4):
        outfh.write(header[i])
    outfh.close()
    geno = new_geno_base+'-credible-wigs-geno.txt'
    os.system('cat %s %s > %s' % (geno_head,geno_body,geno))

def add_pop_to_smartpca(smartpca_inbase,smartpca_outbase,db_name='FULL',prefix='RB'):
    td = get_table_as_dict(db_name,suppress_fc_check=True)
    pop_lookup = dict([(prefix+d['id'],d['location']) for d in td])
    infile = smartpca_inbase+'.ind'
    outfile = smartpca_outbase+'.ind'
    outfh = open(outfile,'w')
    for l in open(infile):
        ind = l.strip().split()[0]
        newl = '\t'.join(l.split()[:2]+[pop_lookup[ind]+'\n'])
        outfh.write(newl)
    outfh.close()

    open(smartpca_outbase+'.par','w').write('genotypename:\t%s.ancestrymapgeno\n' \
                                   'snpname:\t%s.snp\n' \
                                   'indivname:\t%s.ind\n' \
                                   'evecoutname:\t%s.evec\n' \
                                   'evaloutname:\t%s.eval\n' \
                                   'snpweightoutname:\t%s.snpweightout\n'
                                   'phylipoutname:\t%s.fst\n' % tuple([smartpca_outbase]*7))
        
    ret = os.system('cp %s %s' % (smartpca_inbase+'.ancestrymapgeno',smartpca_outbase+'.ancestrymapgeno'))
    if ret != 0:
        raise OSError
    ret = os.system('cp %s %s' % (smartpca_inbase+'.snp',smartpca_outbase+'.snp'))
    if ret != 0:
        raise OSError
    
