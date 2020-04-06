#!/usr/bin/env python
"""
Unwrap fasta file so that each sequence takes up only one line.

Usage:
	%program  input_file  output_file
"""

import os, sys


def read_wrapped_or_unwrapped_fasta (infile):
	
	outlines = []
	with open(infile, "r") as INFILE:
	
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
	
	i=0
	j=1
	out_dict = {}
	for x in range(len(outlines)/2):
		out_dict[outlines[i]] = outlines[j]
		i += 2
		j += 2
	
	return out_dict


def find_input_fastas (pattern):
	
	all_files = os.listdir("./")	
	fastas = [x for x in all_files if pattern in x]
	return fastas


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

pattern = "specific_genome_sequence.fasta"
outfilename = "my.phy.txt"


all_fastas = find_input_fastas(pattern)
fasta_list = {}
all_loci = []
for i in all_fastas:
	fasta_entries = read_wrapped_or_unwrapped_fasta (i)
	fasta_list[i] = fasta_entries 
	all_loci += fasta_entries.keys() 


all_loci = list( set( all_loci ) )
print "loci / alignment partitions:	", len(all_loci)

out_dict = {}
for i in fasta_list.keys():
	concat_seq = ""
	for l in all_loci:
		try:
			concat_seq += fasta_list[i][l]
		except KeyError:
			concat_seq += "-"*int( l.split("L")[1] ) ### luckily we have the loci-lengths encoded in the fasta headers! So we jsut take it from there.
	out_dict[i] = concat_seq

write_phylip (out_dict, outfilename)

