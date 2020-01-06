
import sys

"""
takes a phylip alignment of codons (!) and extracts 4-fold degenerate sites, returns as phylip format!

the aim:

Estimation of Evolutionary Distance.

as in:
Kumar, S. and Subramanian, S. (2002). Mutation rates in mammalian genomes. PNAS 99, 803-808.

methods section:
Evolutionary divergence (d4) between sequences at fourfold-degenerate sites was estimated by using the Tamura-Nei method (22) to correct for multiple hits by accounting for transition/transversion rate and base-frequency biases. 

this "standardised" measure of evolutionary distance will allow to place the Nepenthes radiation into comparison with familiar evolutionary splits (e.g. human - chimp) and the divergence in other radiations (e.g. mammals)

## but maybe not necessary? Here from chimp genome paper:

The Chimpanzee Sequencing and Analysis Consortium (2005). Initial sequence of the chimpanzee genome and comparison with the human genome. Nature 437, 69-87.

comparing ALL sites / whole genomes, coding or not:
"Single-nucleotide substitutions occur at a mean rate of 1.23% between copies of the human and chimpanzee genome, with 1.06% or less corresponding to fixed divergence between the species."

"""

fourfold_degen_codons = ["GCN", "CGN", "GGN", "CTN", "CCN", "TCN", "ACN", "GTN"]
# source: https://en.wikipedia.org/wiki/Codon_degeneracy
# there are 8 4-fold degenerate codons: codons that may contain any base at 1 positions (pos3) and still encode the same amino acid

fourfold_degen_codon_starts = set( [ x.rstrip("N") for x in fourfold_degen_codons ] )

def check_if_fourfold_degenerate (codon_list):
	
	fourfold_degs = [ x for x in codon_list if x[:2] in fourfold_degen_codon_starts ]
#	print fourfold_degs
	if len( fourfold_degs ) == len( [ x for x in codon_list if x != "---"] ): # allows missing codons
		fourfold_deg_column = [ x[2] for x in codon_list ]
#		print "4-fold degenerate"
		return fourfold_deg_column
	else:
#		print "not 4-fold degenerate"
		return None

def read_phylip(INFILE):
	
	indict = {}
	with open(INFILE, "r") as infile:
		infile.readline() # remove header
		for line in infile:
			fields = line.strip("\n").split()
			indict[fields[0]] = fields[1]
	
	for k,v in indict.items():
		print k, len(v)
	
	return indict


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


def walk_alignment (indict):
	
	n_columns = len( indict.values()[0] )
	
	tax_order = indict.keys()
	
	fourfold_columns = []
	idx=0
	while idx < len( indict.values()[0] ):
		codon_column = [ indict[x][idx:(idx+3)] for x in tax_order ]
		fourfold_col = check_if_fourfold_degenerate (codon_column)
		fourfold_columns.append( fourfold_col )
#		print fourfold_col
		idx += 3
	
	good_columns = [ x for x in fourfold_columns if x ]
	print len( good_columns )
	
	out_dict = {tax:[] for tax in tax_order}	
	for idx, tax in enumerate( tax_order ):
		for col in good_columns:
			out_dict[tax].append( col[idx] )
	
	
	out_dict_fourfold_degen_sites = {k: "".join(v) for k,v in out_dict.items() }	
	return out_dict_fourfold_degen_sites

def count_pairwise_naive_distance ( align_dict ):
	
	all_pairs = set()
	for p in align_dict.keys():
		for q in align_dict.keys():
			if p != q:
				all_pairs.add( "$$$".join( sorted([p,q]) ) )
	
	outl = []
	outl.append("/t".join(["tax1","tax2","SNPs_at_4fold_degen_sites","total_4fold_degen_sites","percent_difference_at_4fold_degen_sites"]))
	pairs_dist = {}
	for pair in all_pairs:
		diff = 0
		total = 0
		tax1 = pair.split("$$$")[0]
		tax2 = pair.split("$$$")[1]
		for idx in range(len( align_dict.values()[0] )):
			p = align_dict[tax1][idx]
			q = align_dict[tax2][idx]
			if p != "-" and q != "-":
				total += 1
				if p != q:
					diff += 1
		pairs_dist[pair] = [diff, total, round(float(diff)/float(total)*100, 2) ]
		print pair.split("$$$"), pairs_dist[pair]
		outl.append( "/t".join([tax1,tax2,str(diff),str(total),str(pairs_dist[pair][2]) ]) )
	
	return "\n".join(outl) + "\n"

def count_pairwise_distance_all_and_fourfold_degen ( raw_dict, fourfold_degen_dict ):
	
	all_pairs = set()
	for p in raw_dict.keys():
		for q in raw_dict.keys():
			if p != q:
				all_pairs.add( "$$$".join( sorted([p,q]) ) )
	
	outl = []
	outl.append("\t".join(["tax1","tax2","SNPs_all","total_sites","percent_difference_all_sites","SNPs_at_4fold_degen_sites","total_4fold_degen_sites","percent_difference_at_4fold_degen_sites"]))
	for pair in all_pairs:
		outline = []
		tax1 = pair.split("$$$")[0]
		tax2 = pair.split("$$$")[1]
		outline.append(tax1)
		outline.append(tax2)
		print tax1, tax2
		diff = 0
		total = 0
		print len( raw_dict.values()[0] )
		for idx in range(len( raw_dict.values()[0] )):
			p = raw_dict[tax1][idx]
			q = raw_dict[tax2][idx]
			if p != "-" and q != "-":
				total += 1
				if p != q:
					diff += 1
		outline.append(str(diff))
		outline.append(str(total))
		outline.append(str( round(float(diff)/float(total)*100, 2) ))

		diff = 0
		total = 0
		for idx in range(len( fourfold_degen_dict.values()[0] )):
			p = fourfold_degen_dict[tax1][idx]
			q = fourfold_degen_dict[tax2][idx]
			if p != "-" and q != "-":
				total += 1
				if p != q:
					diff += 1
		outline.append(str(diff))
		outline.append(str(total))
		outline.append(str( round(float(diff)/float(total)*100, 2) ))

		outl.append( "\t".join( outline ) )
		print "\t".join( outline )
		
	return "\n".join(outl) + "\n"
							
	
	
###############################
# main

infilename = sys.argv[1]

in_alignment = read_phylip(infilename)
#print in_alignment
out_dict_fourfold_degen_sites = walk_alignment( in_alignment )
#report = count_pairwise_naive_distance ( in_alignment )
#report = count_pairwise_naive_distance ( out_dict_fourfold_degen_sites )

report = count_pairwise_distance_all_and_fourfold_degen (in_alignment, out_dict_fourfold_degen_sites)

with open(infilename + ".divergence_report.txt", "w") as F:
	F.write(report)


write_phylip (out_dict_fourfold_degen_sites, infilename + ".fourfold_deg_sites.txt")

print "Done!"


