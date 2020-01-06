#!/usr/bin/env python
"""
reads a fasta file and exports a phylip file, nothing else.
"""

import os, sys


def read_wrapped_or_unwrapped_fasta (infile):
	
	outlines = []
	with open(infile, "r") as INFILE:
	
		first_id = INFILE.readline().strip("\n").strip(">")
		outlines.append(first_id)
		seq = ""
	
		for line in INFILE:
			line = line.strip("\n")
			if line.startswith(">"):
				outlines.append(seq)
				outlines.append(line.strip(">").strip("\n"))
				seq = ""
			else:
				if len(line) > 0:
					seq += line.strip("\n")
	
		# append last seq
		outlines.append(seq)
	
	i=0
	j=1
	out_dict = {}
	for x in range(len(outlines)/2):
		out_dict[outlines[i]] = outlines[j]
		i += 2
		j += 2
	
	return out_dict


def write_phylip (concat_dict, outfile):
	
	with open(outfile, "w") as OUTFILE:
		ntaxa = len(concat_dict.keys())
		len_align = len(concat_dict[concat_dict.keys()[0]]) # gets value of first element in dictionary -> this OK since all seqs have same length
		header = str(ntaxa) + " " + str(len_align) + "\n"
		OUTFILE.write(header)
		for sample, seq in concat_dict.items():
			out_line = sample + "    " + seq + "\n"
			OUTFILE.write(out_line)
	OUTFILE.close()

######################

infilename = sys.argv[1]
outfilename = sys.argv[2]

fasta_data = read_wrapped_or_unwrapped_fasta (infilename)

write_phylip (fasta_data, outfilename)

