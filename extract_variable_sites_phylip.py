#!/usr/local/bin/python
# Python 2.7.6 / 2.6
# extract_variable_sites_phylip.py
# 20 Feb 2015
# Mathias Scharmann

"""
extracts the variable positions from a phylip alignment

usage:

toolbox=/gdc_home3/schamath/tools
python $toolbox/extract_variable_sites_phylip.py raff_complex_Brunei_r0.75.phy

"""

# prob. best in python 2.6 since Bio is nor installed for 2.7


from Bio import AlignIO
import sys


def f9(seq):
    # Not order preserving
    return {}.fromkeys(seq).keys()

# get varpos:
def find_varpos_msa_dict(msa_dict):
	the_taxorder = msa_dict.keys()
	the_list = []
	for i in range(0, len(msa_dict[msa_dict.keys()[0]])):
		print i
		column = "".join( [msa_dict[tax][i] for tax in the_taxorder] )
		if len( f9(column.replace("n", "")) ) > 1:
			the_list.append(column)
	out_dict = {}
	for idx, tax in enumerate(the_taxorder):
		print tax
		outstr = ""
		for i in range(0, len(the_list) ): 
			outstr += the_list[i][idx]
		out_dict[tax] = outstr
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

## main

infile=sys.argv[1]
input_handle = open(infile, "rU")

inalign = AlignIO.read(input_handle, "phylip-relaxed")

# get simple dictionary out of the msa object:

msa_dict = {}
for entry in inalign:
	msa_dict[entry.name] = str(entry.seq)

out_dict = find_varpos_msa_dict(msa_dict)

write_phylip (out_dict, "{0}.SNPs.phy".format(infile))
