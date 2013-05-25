#!/usr/bin/env python

from short_read_analysis import preprocess_radtag_lane,extract_genotypes_from_mclgr
import os,sys

def add_pheno_to_map(phenotypes,maploci,genotypes):
    phenomap = {}

    phenoset = set(reduce(lambda x,y:x+y, [rec.keys() for rec in phenotypes]))
    phenoset.discard('id')
    phenomaploci = {}.fromkeys(list(phenoset),('',''))
    phenomaploci.update(maploci)
    
    for pd in phenotypes:
        if pd.has_key('id'):
            if genotypes.has_key(pd['id']):
                #print >> sys.stderr, pd['id']
                pmd = dict([(k,v) for k,v in pd.items() if k != 'id'])
                pmd.update(genotypes[pd['id']])
                phenomap[pd['id']] = pmd
            else:
                print >> sys.stderr, 'no matching genotypes for pheno line %s' % pd['id']
        else:
            print >> sys.stderr, 'no id in %s' % pd
            
    return phenomaploci,phenomap

if __name__ == '__main__':
    db,mapfile,outfile = sys.argv[1:4]

    if ',' in mapfile:
        mapf,mIDf = m.split(',')
    else:
        mapf = mapfile
        mIDf = False

    if ',' in db:
        phenotypes = []
        for db_i in db.split(','):
            phenotypes.extend(preprocess_radtag_lane.get_table_as_dict(db_i,suppress_fc_check=True))
    else:
        phenotypes = preprocess_radtag_lane.get_table_as_dict(db,suppress_fc_check=True)
    
    maploci,genotypes = extract_genotypes_from_mclgr.load_cross_radtag_genotypes(mapf,mIDf)
    
    phenomaploci,phenomap = add_pheno_to_map(phenotypes,maploci,genotypes)
    print >> sys.stderr, '%s pheno+map loci, %s lines' % (len(phenomaploci),len(phenomap))
    og,mID = extract_genotypes_from_mclgr.output_cross_radtag_genotypes(phenomaploci,phenomap,outfile)

    
