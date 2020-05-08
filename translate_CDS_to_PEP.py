
import sys

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
		   "TGG" : "W", "CGG" : "R", "AGG" : "R", "GGG" : "G" 
		   }


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


def translate ( cds_file ):
	
	dna_seqs = read_fasta( cds_file )
	
	out_lines = []
	for id, seq in dna_seqs.items():
		protein_sequence = ""
		for i in range(0, len(seq)-(3+len(seq)%3), 3):
			protein_sequence += codon_table[seq[i:i+3]]
		out_lines.append( ">" + id )
		out_lines.append( protein_sequence )
	
	outname = cds_file + ".translated.fa"
	with open( outname, "w" ) as OUTF:
		OUTF.write(	"\n".join( out_lines ) + "\n")
	
####

if len(sys.argv) != 2:
	print "usage: python translate_CDS_to_PEP.py cds_fasta"
	exit()
	
translate(sys.argv[1])
