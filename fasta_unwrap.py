#!/usr/bin/env python
"""
Unwrap fasta file so that each sequence takes up only one line.

Usage:
	%program  input_file  output_file
"""

import sys

# try:
# 	from Bio import SeqIO
# except:
# 	print "This program requires the Biopython library"
# 	sys.exit(0)
# 
# try:
# 	in_file = open(sys.argv[1], "rU")
# 	out_file = open(sys.argv[2], "w")
# except:
# 	print __doc__
# 	sys.exit(0)
# 
# sequences = ([seq.id, seq.seq.tostring()] for seq in SeqIO.parse(in_file, 'fasta'))
# with open(sys.argv[2], "w") as out_file:
# 	for seq in sequences:
# 		out_file.write(">" + seq[0] + "\n" + seq[1] + "\n")



#####


outlines = []

with open(sys.argv[1], "r") as INFILE:
	
	first_id = INFILE.readline().strip("\n")
	outlines.append(first_id)
	seq = ""
	
	for line in INFILE:
		line_clean = line.strip("\n")
		if line_clean.startswith(">"):
			outlines.append(seq)
			outlines.append(line_clean)
			seq = ""
		else:
			if len(line_clean) > 0:
				seq += line_clean
	
	# append last seq
	outlines.append(seq)
	
with open(sys.argv[2], "w") as OUTFILE:
	OUTFILE.write("\n".join(outlines))
	