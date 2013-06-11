'''functionality for enclosure experiments

prefix refers to string identifying indivs in sequencing samples (and therefore vcf files)

includes:
generate phenotype files and handle generation of genotypes (using variant_detection.write_wigs_output)
'''
try:
    from rtd.preprocess_radtag_lane import get_table_as_dict
except:
    from radtag_denovo.preprocess_radtag_lane import get_table_as_dict

def indivs_in_vcf(vcf_obj):
    indivs = set([])
    for v in vcf_obj.values():
        these_indivs = set(v['indiv_gt'].keys())
        indivs = indivs.union(these_indivs)
    return indivs

def indiv_lists_by_enc(vcf_obj,db_name,enc,sites=('D','L'),prefix='RB',sort_key=None):
    td = get_table_as_dict(db_name,suppress_fc_check=True,sq='enc=%s' % enc)
    vcf_indivs = indivs_in_vcf(vcf_obj)
    indiv_lists = tuple([[d['id'] for d in td if d['site'] == this_site and prefix+d['id'] in vcf_indivs] for this_site in sites])
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

def write_wigs_all_simple(vcf_obj,db_name,enc,site,outbase):
    indiv_lists = indiv_lists_by_enc(vcf_obj,db_name,enc,sites=[site])
    phenos = wigs_phenotypes_from_DB(db_name,indiv_lists)

    open(outbase+'-wigs-indivs.tuple','w').write(indiv_lists.__repr__())
    write_wigs_pheno(phenos,outbase+'-wigs-pheno.txt')

    write_wigs_genotypes(vcf_obj,indiv_lists,outbase,[0],[0],'RB')

def cut_fn(sd): #NOTE QD 5 DIFFERS FROM SINERGIA/CRL
    summ_stats = ['FS','QUAL','BaseQRankSum','QD','SB','ReadPosRankSum']
    FS,QUAL,BaseQRankSum,QD,SB,ReadPosRankSum = map(float,[sd.get(ss,0) for ss in summ_stats])
    return FS < 72 and QUAL > 218 and BaseQRankSum < 5 and QD >= 5 and SB < 10 and ReadPosRankSum > -9 and sd['fh'] < 0.55
