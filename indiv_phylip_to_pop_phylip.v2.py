#!/usr/local/bin/python
# Python 2.7.6
# 
# 2 Mar 2014
# Mathias Scharmann


import sys

revcombs = {'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'K': 'GT', 'M': 'AC', 'N': 'ACGT', 'S': 'CG', 'R': 'AG', 'W': 'AT', 'V': 'ACG', 'Y': 'CT', 'A' : 'A', 'C':'C', 'T':'T', 'G':'G', 'n':'-', '-':'-'}

def get_IUPAC_amb_code (nucs):
	
	# nucs must be a "".joined string of the sorted list of nucleotides
	
	# the IUPAC ambiguity codes:
	# the keys in this IUPAC ambiguity code dictionary are produced by:
	# "".join(sorted(list("NUCLEOTIDES"))) -> consistent lookup possible without additional loops!
	combs = {'AC':'M', 'GT':'K', 'CG':'S', 'AT':'W', 'AG':'R', 'CT':'Y',
	'ACG':'V', 'CGT':'B', 'AGT':'D', 'ACT':'H', 'ACGT':'N',
	'A' : 'A', 'C':'C', 'T':'T', 'G':'G', 'n':'-', '-':'-', 'N':'N'}

	
	out_code = combs[nucs]
	
	return out_code
	
	
def read_phylip(INFILE, samples):
	
	indict = {}
	with open(INFILE, "r") as infile:
		infile.readline() # remove header
		for line in infile:
			fields = line.strip("\n").split("    ")
			if fields[0] in samples: # discard unwanted samples straight away
				indict[fields[0]] = fields[1]

	return indict

def read_pop_info (INFILE):
	
	popdict = {}
	with open(INFILE, "r") as infile:
		for line in infile:
			fields = line.strip("\n").split("\t")
			if fields[1] in popdict.keys():
				popdict[fields[1]].append(fields[0])
			else: 
				popdict[fields[1]] = [fields[0]]

	print popdict
	return popdict
	
def make_pop_level_seqs (indict, popdict):
	
	# sort samples in indict into their populations:
	binned_dict = {}
	for sample, seq in indict.items():
		for pop, samples in popdict.items():
			if sample in samples:
				try:
					binned_dict[pop].append(seq)
				except KeyError:
					binned_dict[pop] = [seq]
	
	# make consensus sequence out of individual sequences in the binned_dict:
	consensus_dict = {}
	for pop, seqs in binned_dict.items():
		outseq = []
		print "working on pop {0}".format(pop)
		for column_index in range(0, len(seqs[0])): # to loop over all columns
			aligned_column = [ seq[column_index] for seq in seqs ]
			resolved = [ revcombs[x] for x in aligned_column ] 
			aligned_column = list("".join(resolved))
			nucs = [x for x in aligned_column if x != '-'] # remove all the n -
			if len(nucs) == 0:
				outnuc = "-"
			else:
				outnuc = get_IUPAC_amb_code( "".join(sorted(list(set(nucs)))) )
			print column_index, "\r",	
			outseq.append(outnuc)
		outseqstring = "".join(outseq)
		consensus_dict[pop] = outseqstring
	

	return consensus_dict

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

################################## MAIN

def main(argv=None):
	if argv is None:
		argv = sys.argv[1:]
	
	popdict = read_pop_info(argv[1])
	samples = []
	for x in popdict.values():
		samples += x
	print samples
	indict = read_phylip(argv[0], samples)
	consensus_dict = make_pop_level_seqs (indict, popdict)
	write_phylip (consensus_dict, argv[2])
	
	print "Done!"
	
	
if __name__ == "__main__":
	main(sys.argv[1:])
