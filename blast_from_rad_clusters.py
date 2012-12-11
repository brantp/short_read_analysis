#!/usr/bin/env python

import os,sys
from short_read_analysis import extract_genotypes_from_mclgr

def get_cluster_seq(clid,clusters,labels,seqindex=0):
    if '.' in clid:
        clid = clid.split('.')[0]
    return sorted([(int(labels[lid].split('.')[1]),labels[lid].split('.')[2]) for lid in clusters[clid]],reverse=True)[seqindex][1]

def get_top_reads(seq,uniquedlines):
    hits = []
    for l in uniquedlines:
        if l.startswith(seq):
            hits.append(l.strip())

    #return hits
    try:
        r1_pick,n,qual,ind,ct,r2,r2ct = sorted([(int(l.split()[1]), l.split()) for l in hits],reverse=True)[0][1]
        r2_pick = sorted(zip(r2ct.split(','),r2.split(',')),reverse=True)[0][1]
    except:
        return None,None

    return r1_pick,r2_pick

if __name__ == "__main__":

    gr,tab,uniqued,seqfile,blastdb,blast_out = sys.argv[1:7]
    clids = [l.strip() for l in sys.stdin]

    outroot,grname = os.path.split(gr)
    #outroot = '/scratch/brantp'

    if os.path.exists(seqfile):
        print >> sys.stderr, 'sequence file %s present, proceed' % seqfile
    else:
        clusters,labels,mID_by_label = extract_genotypes_from_mclgr.load_cluster_data(gr,tab)
        print >> sys.stderr, 'load uniqued lines from %s ...' % uniqued,
        uniquedlines = open(uniqued).readlines()
        print >> sys.stderr, 'done'
        
        print >> sys.stderr, 'write longest pairs from %s for clusters %s/%s to %s' % (uniqued,gr,tab,seqfile)

        sfh = open(seqfile,'w')

        for clid in clids:
            print >> sys.stderr, clid,

            for i in range(2):
                seq = get_cluster_seq(clid,clusters,labels,seqindex=i)
                r1,r2 = get_top_reads(seq,uniquedlines)

                if r1 is not None and r2 is not None and r2 is not '.':
                    print >> sys.stderr, i,'f/r',
                    sfh.write('>%s-%sf\n%s\n>%s-%sr\n%s\n' % (clid,i,r1,clid,i,r2))

            print >> sys.stderr, ''

        sfh.close()

    
    if os.path.exists(blast_out):
        print >> sys.stderr, 'blast output %s present, proceed' % blast_out
    else:
        print >> sys.stderr, 'run blast %s against %s, write to %s' % (seqfile,blastdb,blast_out)
        os.system('blastall -p blastn -e 0.01 -m 9 -d %s -i %s > %s' % (blastdb,seqfile,blast_out))
    
    print >> sys.stderr, 'done'
            
    

    
