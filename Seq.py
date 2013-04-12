from paths import paths

from UserDict import UserDict
from UserString import UserString
from copy import deepcopy
from glob import glob
import os, sys, re, cPickle, subprocess, shutil, random, Util

gencode = {
    'ATA':'I', #Isoleucine
    'ATC':'I', #Isoleucine
    'ATT':'I', # Isoleucine
    'ATG':'M', # Methionine
    'ACA':'T', # Threonine
    'ACC':'T', # Threonine
    'ACG':'T', # Threonine
    'ACT':'T', # Threonine
    'AAC':'N', # Asparagine
    'AAT':'N', # Asparagine
    'AAA':'K', # Lysine
    'AAG':'K', # Lysine
    'AGC':'S', # Serine
    'AGT':'S', # Serine
    'AGA':'R', # Arginine
    'AGG':'R', # Arginine
    'CTA':'L', # Leucine
    'CTC':'L', # Leucine
    'CTG':'L', # Leucine
    'CTT':'L', # Leucine
    'CCA':'P', # Proline
    'CCC':'P', # Proline
    'CCG':'P', # Proline
    'CCT':'P', # Proline
    'CAC':'H', # Histidine
    'CAT':'H', # Histidine
    'CAA':'Q', # Glutamine
    'CAG':'Q', # Glutamine
    'CGA':'R', # Arginine
    'CGC':'R', # Arginine
    'CGG':'R', # Arginine
    'CGT':'R', # Arginine
    'GTA':'V', # Valine
    'GTC':'V', # Valine
    'GTG':'V', # Valine
    'GTT':'V', # Valine
    'GCA':'A', # Alanine
    'GCC':'A', # Alanine
    'GCG':'A', # Alanine
    'GCT':'A', # Alanine
    'GAC':'D', # Aspartic Acid
    'GAT':'D',    # Aspartic Acid
    'GAA':'E',    # Glutamic Acid
    'GAG':'E',    # Glutamic Acid
    'GGA':'G',    # Glycine
    'GGC':'G',    # Glycine
    'GGG':'G',    # Glycine
    'GGT':'G',    # Glycine
    'TCA':'S',    # Serine
    'TCC':'S',    # Serine
    'TCG':'S',    # Serine
    'TCT':'S',    # Serine
    'TTC':'F',    # Phenylalanine
    'TTT':'F',    # Phenylalanine
    'TTA':'L',    # Leucine
    'TTG':'L',    # Leucine
    'TAC':'Y',    # Tyrosine
    'TAT':'Y',    # Tyrosine
    'TAA':'_',    # Stop
    'TAG':'_',    # Stop
    'TGC':'C',    # Cysteine
    'TGT':'C',    # Cysteine
    'TGA':'_',    # Stop
    'TGG':'W',    # Tryptophan
    }

degen = {'B': ['C', 'G', 'T'],
	 'D': ['A', 'G', 'T'],
	 'H': ['A', 'C', 'T'],
	 'K': ['G', 'T'],
	 'M': ['A', 'C'],
	 'R': ['A', 'G'],
	 'S': ['C', 'G'],
	 'V': ['A', 'C', 'G'],
	 'W': ['A', 'T'],
	 'Y': ['C', 'T']}


class Sequence(UserString):
    "simple string-descended sequence class."

    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N','X': 'X', 'R':'Y', 'Y':'R', 'K':'M', 'M':'K', 'W':'S', 'S':'W','H':'D','V':'B','D':'H','B':'V'}

    def __init__(self, seq):
        UserString.__init__(self,seq)

    def rc(self):
        "returns the reverse complement of stored sequence"
        sl = list(self.data)
        sl.reverse()
        return Sequence(''.join([self.complement[n.upper()] for n in sl]))

    def undegenerate(self,include_n=True,return_both=False):
        '''calls undegenerate_iupac() on all bases in self, choosing one at random for each degenerate base

        setting include_n=False will pass Ns thru unresolved (i.e. masked chars stay masked)
        return_both is a special runmode in which the first two bases for an ambiguity code are returned as a pair of sequences
        (this is largely for writing unphased haplotypes)
        '''

        if return_both:
            s1 = []
            s2 = []
            for c in self.data:
                ui = undegenerate_iupac(c,False)
                if len(ui) == 1:
                    s1.append(ui[0])
                    s2.append(ui[0])
                else:
                    s1.append(ui[0])
                    s2.append(ui[1])
            return ''.join(s1), ''.join(s2)
        else:
            return ''.join([random.choice(undegenerate_iupac(c,include_n)) for c in self.data])
        

    def translate_sequence(self,frame=0):
        '''takes a sequence, and translates it after parsing it into codon triplets
        returns the translation'''

        codonLi = re.findall(r'[actg]{3}',self.data[frame:],re.I)
        aaLi = [gencode[cod] for cod in codonLi]

        return ''.join(aaLi[:-1]) #don't give back the stop!

    def longest_orf(self,forward=True,format='frame'):
        '''returns frame of longest uninterrupted translation
        currently only works on plus strand, returns sequence'''

        trans = [(len(self.translate_sequence(f)),f) for f in range(3)]
        if format == 'frame':
            return max(trans)[1]
        elif format == 'seq':
            seqs_lens = [(len(s),s) for s in self.translate_sequence(max(trans)[1]).split('_')]
            return max(seqs_lens)[1]

class Fasta(UserDict):
    '''Store and operate on sequences in FASTA format

    init argument is filename; if specified loads from fasta (otherwise call load)

    seqs_to_load is a list of strings to be evaluated as regex against sequence names
    if present, only sequences matching at least one regex string will be loaded

    printing or assigning from (e.g. __repr__) produces fasta format'''

    def __init__(self, filename=None, seqs_to_load=None, seqobj=Sequence):
        UserDict.__init__(self)
        self.order = []
        self.metadata = {}
        if filename: 
            self.load(filename,seqs_to_load,seqobj)
        else:
            self.filename = None

    def invalid_chars(self,line,chars):
        '''returns list of invalid chars in line (nonmatches to list chars)'''

        return list(set([c for c in list(line) if c not in chars]))

    def load(self, filename, seqs_to_load, seqobj):
        "load sequences from fasta file"

        if isinstance(filename,str):
            if os.path.exists(filename.rstrip('.gz')+'.pkl') and seqs_to_load is None:
                self.update(cPickle.load(open(filename.rstrip('.gz')+'.pkl')))
                self.filename = filename

            else:
                import gzip
                if filename.endswith('.gz'):
                    opener = gzip.open
                else:
                    opener = open
        
                fsock = opener(filename)
        else:
            fsock = filename
            
        seq_name = None
        seq_read = []
        
        for line in fsock.readlines():
            if line.startswith(">"):
                if seq_name is not None:
                    self[seq_name] = seqobj(''.join(seq_read))
                if seqs_to_load:
                    if not any([re.search(s,line) for s in seqs_to_load]):
                        seq_name = None
                        continue
                seq_name = line.strip("\n")[1:line.find(' ')>0 and line.find(' ') or None]
                if line.find(' ')>0: self.metadata[seq_name] = line.strip("\n")[line.find(' ')+1:]
                seq_read = []
                self.order.append(seq_name)
            elif len(line.strip()) > 0:
                if seqs_to_load:
                    if not any([re.search(s,str(seq_name)) for s in seqs_to_load]):
                        continue
                invchars = self.invalid_chars(line.strip("\n").upper(), seqobj.complement.keys())
                invchars = False
                if invchars:
                    raise ValueError, "not valid: "+','.join(invchars)
                else:
                    seq_read.append(line.strip("\n").upper())
        if seq_name: self[seq_name] = seqobj(''.join(seq_read))

        self.filename = filename
        fsock.close()


    def __repr__(self):
        '''returns sequences in fasta format, e.g.:

        >name
        sequence
        '''
        return "\n".join([">%s\n%s" % (k, v) for k, v in self.data.items() if k[0:2] != "__"])

    def repr_opts(self,ids=None,metadata=1,line_length=None):
        '''since __repr__ blindly gives back the whole Fasta, fasta stylee,

        repr_opts gives some choices: list "ids" returns only sequences in the list
        metadata=1 puts metadata (stuff after a space in header) back

        if line_length is set, adds a newline every line_length characters in sequence output
        '''

        ids = ids or list(set(self.seq_names()))

        if not isinstance(ids, list):
            ids = [ids]

        if metadata and self.metadata: 
            data = []
            for k in ids:
                try:
                    data.append(">%s %s\n%s" % (k, self.metadata[k], self[k]))
                except KeyError:
                    data.append(">%s\n%s" % (k, self[k]))
        else:
            data = [">%s\n%s" % (k, self[k]) for k in ids]
        if line_length is not None:
            lbdata = []
            for d in data:
                h,s = d.split('\n')
                lis = list(s)
                for i in range(0, len(s)+(len(s)/line_length)+1,line_length):
                    lis[i:i] = '\n'
                lbdata.append(h+'\n'+''.join(lis))
            data = lbdata
        return '%s\n' % "\n".join(data)

    def seq_len(self,include_seq=None):
        '''return lengths of each sequence as dict.

        set list "include_seq" to return only some seqs'''
        include_seq = include_seq and include_seq or self.seq_names()
        return dict([(k, len(v)) for k, v in self.data.items() if k in include_seq])

    def rc_sequences(self,ids=None):
        '''revcomps sequences in list ids. 

        defaults to all seqs'''
        ids = ids and ids or self.seq_names() 
        for id in ids:
            self[id] = self[id].rc()

    def longest_orfs(self):
		'''returns new Seq.Fasta object of longest orfs in each seq'''
		new = Fasta()
		new.metadata = deepcopy(self.metadata)
		for k,v in self.items():
			new[k] = v.longest_orf(format='seq')
		return new

    def write_to_file(self, outfile=None, set_filename=0, ids=None, metadata=1, line_length=None):
        '''output current fasta to specified file; if no file is specified, self.filename is used 

        if set_filename is specified, changes the filename of this instance

        defaults to all sequences, if list ids specified, only output a subset'''

        outfile = (outfile or self.filename)

        try:
            os.unlink(outfile)
        except OSError:
            pass

        fsock = open(outfile, "w")
        fsock.write(self.repr_opts(ids,metadata,line_length))
        fsock.close()

        if set_filename:
            self.filename = outfile

    def seq_names(self):
        '''returns list of all sequence names in Fasta object

        '''
        if self.order:
            return [k for k in self.order if k in self.data.keys()] + list(set(self.data.keys())-set(self.order))
        else:
            return [k for k in self.keys() if k[0:2] != "__"]


    def split_multifasta(self,directory=None,basenames_only=1,ext='.fa'):
        '''if directory specified, deposits single sequence fastas in directory
        if basenames_only=1, return only basenames, otherwise return entire path to files

        if no directory is supplied:
        returns a list of Fasta objects containing a single sequence each <- NOT YET IMPLEMENTED'''

        if directory:
            files = []
            os.path.exists(directory) or os.makedirs(directory)
            for seq in self.seq_names():
                filename = seq+ext
                outfile = os.path.join(directory,filename)
                self.write_to_file(outfile,ids=[seq])
                if basenames_only:
                    files.append(filename)
                else:
                    files.append(outfile)
            return files
                

    def substr_from_gff(self, region_list, name_key="Name", use_seqid=True, plus_strand=False):
        '''given a list of gff regions, return a new Fasta object.

        sequence names are GFF.Region["seqid"]+'.'+GFF.Region["attributes"][name_key] (if use_seqid=1)

        at present, duplicate final labels will be overwritten

        if plus_strand=1, sequence is returned on plus strand ( i.e. seq from '-' features is rc()\'d )'''
        if not isinstance(region_list, list):
            region_list = [region_list]

        new_fasta = Fasta()
        for region in region_list:
            if self[region["seqid"]]:
                label = (use_seqid and region["seqid"] or '') + \
                    ((use_seqid and name_key) and '.' or '') + \
                    (name_key and str(region["attributes"][name_key]) or '')
                firstbase = region["start"] - 1
                new_fasta[label] = self[region["seqid"]][firstbase:region["end"]]
                if plus_strand and region["strand"] == '-': 
                    new_fasta[label] = new_fasta[label].rc()
        
        return new_fasta

#FASTA HELPERS
def write_fastas_to_files(fastas,clobber=False):
	'''writes fastas if they\'re instances
	if clobber=False, only writes if .filename doesn\'t exist
	returns list of filenames
	'''
	filenames = []
	for fasta in fastas:
		if isinstance(fasta,str):
			filenames.append(fasta)
		else:
			if clobber or not os.path.exists(fasta.filename):
				fasta.write_to_file()
			filenames.append(fasta.filename)
	return filenames

def num_seqs_in_fasta_file(filename):
    if filename.endswith('.gz'):
        return int(subprocess.Popen('gunzip -c %s | grep ">" | wc' % filename,shell=True,stdout=subprocess.PIPE).stdout.read().split()[0])
    else:
        return int(subprocess.Popen('cat %s | grep ">" | wc' % filename,shell=True,stdout=subprocess.PIPE).stdout.read().split()[0])

def seqnames_from_fasta_file(filename):
    if filename.endswith('.gz'):
        return [l.split()[0].lstrip('>') for l in subprocess.Popen('gunzip -c %s | grep ">" ' % filename,shell=True,stdout=subprocess.PIPE).stdout.readlines()]
    else:
        return [l.split()[0].lstrip('>') for l in subprocess.Popen('cat %s | grep ">" ' % filename,shell=True,stdout=subprocess.PIPE).stdout.readlines()]
    

def split_multifasta_to_dir(fasta,clobber=False):
    '''takes a filename, creates a directory (minus extension, if .fa[sta]) of single files named after seqids

    returns paths to created files'''
    match = re.search(r'(.+?)\.f[nsta]{1,4}(\..+?|$)',fasta)

    if match:
        tdir = match.groups()[0]
    else:
        tdir = fasta

    if clobber:
        shutil.rmtree(tdir)
    
    if os.path.isdir(tdir) and num_seqs_in_fasta_file(fasta) == len(glob(tdir+'/*.fa')):
        return glob(tdir+'/*.fa')
    else:
        os.makedirs(tdir) #if this fails, dir exists
        fa_inst = Fasta(fasta)
        return fa_inst.split_multifasta(tdir,basenames_only=False)
    

def undegenerate_iupac(s,include_n=True):
	'''returns a list of the explicit DNA bases (actg) from an iupac degenerate symbol'''
	
        thisdegen = degen.copy()
        if include_n:
            thisdegen['N'] = ['A', 'C', 'G', 'T']
            
	return thisdegen.get(s.upper(),[s.upper()])

def degenerate_iupac(lol):
    '''given a list of lists where each inner list is all bases at a given position (and outer list indices are positions, ie aln columns)
    returns an iupac consensus'''

    thisdegen = degen.copy()
    thisdegen['N'] = ['A', 'C', 'G', 'T']
    degdict = dict([(tuple(sorted(v)),k) for k,v in thisdegen.items()])
    degdict[('A',)] = 'A'
    degdict[('C',)] = 'C'
    degdict[('G',)] = 'G'
    degdict[('T',)] = 'T'

    consli = [degdict[tuple(sorted(list(set([b.upper() for b in bases]))))] for bases in lol]
    return ''.join(consli)

def peptide_to_degenerate_iupac(pep):
    '''given a peptide sequence returns the iupac consesnsus in DNA for all possible coding sequences'''

    ctab = Util.invert_dict(gencode)
    
    return degenerate_iupac(reduce(lambda x,y:[xel for xel in x]+[yel for yel in y],[numpy.array([list(codon) for codon in codons]).transpose() for codons in [ctab[c] for c in pep]]))

def is_simple(s,size=[1,2,3]):
	'''returns True if the string s is a simple repeat of any unit size in <size>'''
	
	for w in sorted(size):
		if len(set([s[i:i+w] for i in range(0,len(s) - len(s) % w,w)])) == 1:
			return True

class Eland():
	
	def __init__(self, filename=None):
		self.seqs = {}
		self.quals = {}
		self.filename = None
		if filename:
			self.load(filename)
	
	def load(self,filename):
		self.filename = filename
		next = None
		key = None
		for l in open(filename).readlines():
			if l.startswith('@'):
				key = l[1:].strip()
				next = 'seq'
			elif l.startswith('+'):
				key = l[1:].strip()
				next = 'qual'
			elif next == 'seq':
				self.seqs[key] = l.strip()
			elif next == 'qual':
				self.quals[key] = eland_qual_to_int(l.strip())
				
import numpy					
def cumulative_quality(qual):
	return [numpy.mean(qual[:i+1]) for i,v in enumerate(qual)]


def quality_cutoff(qual,cut):
	belowcut = list(numpy.array(cumulative_quality(qual))<cut)
	if any(belowcut):
		return belowcut.index(True)
	else:
		return len(qual)
	
def eland_qual_to_int(qualstr):
	'''given a quality string in ASCII+64, returns a list of integer quality scores'''
	return [ord(i)-64 for i in qualstr]
	
if __name__ == "__main__":

    #test Fasta class (loading, editing, writing)
    test_fasta = Fasta("/home/brant/py_util/unit_test_data/seq.fasta")

    print "%s\n%s" % (test_fasta.filename, test_fasta.seq_len())
    for k in test_fasta.iterkeys():
        test_fasta[k] += "TGGCG"
    test_fasta.write_to_file("/home/brant/temp/temp.fa", 1)    

    print "%s\n%s" % (test_fasta.filename, test_fasta.seq_len())

    other_test_fasta = Fasta("/home/brant/temp/temp.fa")
    print other_test_fasta.seq_len()
    #end Fasta test
    print other_test_fasta.order

    print "test substr_from_gff\n"
    import GFF                
    seqfile = os.path.join(paths['py_testdata'],"eve.ceratitis_capitata.fa")
    gfffile = os.path.join(paths['py_testdata'],"eve.ceratitis_capitata.fa.gff3")
    seq = Fasta(seqfile)
    gff = GFF.File(gfffile)
    evegene = seq.substr_from_gff([region for region in gff
                                  if 'gene_name' in region['attributes'].keys() and region['attributes']['gene_name']=='eve'],
                                  name_key='gene_name',plus_strand=1)
    print evegene
    
