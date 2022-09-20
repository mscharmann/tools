"""
- reads a short read assembly in fasta format and concatenates contigs in order of appearance up to a given maximum size.
-> creates "mock chromosomes", from hundreds of thousands of tiny pieces to a few large pieces; this greatly facilitates read mapping, variant calling, plotting etc.
- gaps in between contigs are filled with 500 "N"
"""

import sys

def read_wrapped_or_unwrapped_fasta (infile):
	
	outlines = []
	seq = ""
	with open(infile, "r") as INFILE:
		for line in INFILE:
			line = line.strip("\n")
			if line.startswith(">"):
				outlines.append(seq)
				seq = ""
			else:
				if len(line) > 0:
					seq += line.strip("\n")
		# append last seq
		outlines.append(seq)
	
	return outlines[1:]




infile = sys.argv[1]
mock_chrom_size = int(sys.argv[2])

inseqs = read_wrapped_or_unwrapped_fasta (infile)

outfastalines = []
outseq = ""
chrom_count = 1
for s in inseqs:
	if not len(outseq) > mock_chrom_size:
		outseq += s
		outseq += "N"*500
	else:
		outfastalines.append(">mock_chrom_" + str(chrom_count))
		chrom_count += 1
		print(chrom_count)
		outfastalines.append(outseq)
		outseq = ""
		

outfastalines.append(">mock_chrom_" + str(chrom_count))
chrom_count += 1
outfastalines.append(outseq)
outseq = ""


with open(infile + ".mock_chroms.fa", "w") as O:
	O.write("\n".join(outfastalines)+"\n")
	

