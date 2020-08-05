#!/usr/local/bin/python

"""
- drops taxa that lost almost everyting

"""

import sys


def read_phylip(INFILE):
	
	indict = {}
	with open(INFILE, "r") as infile:
		infile.readline() # remove header
		for line in infile:
			if len(line) > 1:
				fields = line.strip("\n").split()
				indict[fields[0]] = fields[1]

	return indict


def read_wrapped_or_unwrapped_fasta (infile):
	
	outlines = []
	with open(infile, "r") as INFILE:
		
		first_id = INFILE.readline().strip("\n")
		outlines.append(first_id.strip(">"))
		seq = ""
	
		for line in INFILE:
			line_clean = line.strip("\n")
			if line_clean.startswith(">"):
				outlines.append(seq)
				outlines.append(line_clean.strip(">"))
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


# get varpos:
def find_varpos_msa_dict(msa_dict):
	the_taxorder = msa_dict.keys()
	the_list = []
	for i in range(0, len(msa_dict[msa_dict.keys()[0]])):
		column = "".join( [msa_dict[tax][i] for tax in the_taxorder] )
		if len( set(column.replace("-", "")) ) > 1:
			the_list.append(column)
	out_dict = {}
	for idx, tax in enumerate(the_taxorder):
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


def write_fasta (msa_dict, outfilename):
	
	outlines = []
	for id in sorted(msa_dict.keys()):
		outlines.append( ">"+id+"\n"+msa_dict[id] )
		
	with open(outfilename, "w") as OUTFILE:
		OUTFILE.write("\n".join(outlines) + "\n")


## main

file_format = sys.argv[1]

if not file_format in ["fasta", "phylip"]:
	print "allowed file formats are 'fasta' and 'phylip'"
	exit()

infile=sys.argv[2]

if file_format == "phylip":
	input = read_phylip(infile)
else:
	input = read_wrapped_or_unwrapped_fasta(infile)


out_dict = find_varpos_msa_dict(input)

print "input sites:	", len(input.values()[0])
print "output variable sites:	", len(out_dict.values()[0])

if file_format == "phylip":
	write_phylip (out_dict, "{0}.varpos.phy".format(infile))
else:
	write_fasta (out_dict, "{0}.varpos.fasta".format(infile))

print "Done!"
