
import sys, os

codon_table = {"TTT" : "F", "CTT" : "L", "ATT" : "I", "GTT" : "V",
		   "TTC" : "F", "CTC" : "L", "ATC" : "I", "GTC" : "V",
		   "TTA" : "L", "CTA" : "L", "ATA" : "I", "GTA" : "V",
		   "TTG" : "L", "CTG" : "L", "ATG" : "M", "GTG" : "V",
		   "TCT" : "S", "CCT" : "P", "ACT" : "T", "GCT" : "A",
		   "TCC" : "S", "CCC" : "P", "ACC" : "T", "GCC" : "A",
		   "TCA" : "S", "CCA" : "P", "ACA" : "T", "GCA" : "A",
		   "TCG" : "S", "CCG" : "P", "ACG" : "T", "GCG" : "A",
		   "TAT" : "Y", "CAT" : "H", "AAT" : "N", "GAT" : "D",
		   "TAC" : "Y", "CAC" : "H", "AAC" : "N", "GAC" : "D",
		   "TAA" : "*", "CAA" : "Q", "AAA" : "K", "GAA" : "E",
		   "TAG" : "*", "CAG" : "Q", "AAG" : "K", "GAG" : "E",
		   "TGT" : "C", "CGT" : "R", "AGT" : "S", "GGT" : "G",
		   "TGC" : "C", "CGC" : "R", "AGC" : "S", "GGC" : "G",
		   "TGA" : "*", "CGA" : "R", "AGA" : "R", "GGA" : "G",
		   "TGG" : "W", "CGG" : "R", "AGG" : "R", "GGG" : "G",
		   "---" : "-" }


def read_fasta (infile):

	fasta_dict = {}
	seq = ""
	name = "dummy"
	with open(infile, "r") as F:
		for line in F:
			if line.startswith(">"):
				fasta_dict[name] = seq
				name = line.lstrip(">").rstrip("\n")
				seq = ""
			else:
				seq += line.strip("\n")

	# last record:
	fasta_dict[name] = seq
	del fasta_dict["dummy"]
	return fasta_dict


def parse_and_translate ( cds_file ):
	
	dna_seqs = read_fasta( cds_file )
	
	all_seqs = dna_seqs.values()
	
	prot_seqs = []
	for seq in dna_seqs.values():
		protein_sequence = ""
		for i in range(0, len(seq)-(3+len(seq)%3), 3):
			protein_sequence += codon_table[seq[i:i+3]]
		prot_seqs.append( protein_sequence )
	
	return prot_seqs


def make_maj_consensus (protseqs):
	
	consensus_seq = ""
	for i in range(len(protseqs[0])):
		AAs = [x[i] for x in protseqs if not x[i] == "-"]
		unique_AAs = list(set(AAs))
		counts = []
		for A in unique_AAs:
			counts.append( AAs.count(A) )
		maxidx = counts.index(max(counts))
		consensus_seq += unique_AAs[maxidx]
	return consensus_seq 
	
####

if len(sys.argv) != 2:
	print ( "usage: python translate_CDS_to_PEP.py cds_fasta_dir" )
	exit()

infiles = [x for x in os.listdir(sys.argv[1]) if x. endswith(".CDS.aln")]

outlines = []
cnt = 0
for f in infiles:
	cnt += 1
	print(cnt)
	PEPseqs = parse_and_translate(sys.argv[1] + "/" + f)
	if len(PEPseqs) > 0:
		cons = make_maj_consensus(PEPseqs)
		outlines.append(">" + f.strip(".CDS.aln"))
		outlines.append(cons)
	else:
		with open("maj_cons_problems_log.txt", "a") as O:
			O.write(f + "\n")
			print ( f )
	
with open("maj_cons.fasta", "w") as O:
	O.write("\n".join(outlines) + "\n")

	
	


