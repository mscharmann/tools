## filter_longest_isoform.py

import sys

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
	for x in range(int(len(outlines)/2.0)):
		out_dict[outlines[i]] = outlines[j]
		i += 2
		j += 2
	
	return out_dict



#i = "candidates.pep"

indict = read_wrapped_or_unwrapped_fasta(sys.argv[1])
# indict = read_wrapped_or_unwrapped_fasta(i)


inseqs_ordered = list(indict.keys())
genes = [ "".join(x.split(".")[:-1]) for x in inseqs_ordered]

isoform_dict = {}
for idx,gene in enumerate(genes):
	try:
		isoform_dict[gene].append(inseqs_ordered[idx])
	except KeyError:
		isoform_dict[gene] = [ inseqs_ordered[idx] ]
		
	

outlines = []
for k,v in isoform_dict.items():
	isolengths = [len(indict[i]) for i in v]
	longest_iso = v[ isolengths.index(max(isolengths)) ]
	print (isolengths, v, longest_iso)
	outlines.append(">" + k)
	outlines.append(indict[longest_iso])
	
	
with open(sys.argv[1] + ".longest_isoform.fa", "w") as O:
	O.write("\n".join(outlines)+"\n")
	
	
	