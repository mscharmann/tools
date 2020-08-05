"""
python 2.7

Mathias Scharmann
2020-06-02

counts the numbers of potentially synonymous and potentially nonsynonymous sites (Nei and Gojobori 1986)

Nei, M., Gojobori, T., 1986. Simple methods for estimating the numbers of synonymous and nonsynonymous nucleotide substitutions. Mol. Biol. Evol. 3, 418-426. https://doi.org/10.1093/oxfordjournals.molbev.a040410

it is a 'mutational opportunity method' (Bierne and Eyre-Walker, 2003, Genetics 165: 1587-1597)

the relevant paragraph in Nei & Gojobori (1986) is:

We now compute the number of synonymous sites (s) and the number of non-synonymous sites (n) for each codon, considering the above property of codon changes. We denote by f(i) the fraction of synonymous changes at the ith position of a given codon(i=1,2,3). The s and n for this codon are then given by s = sum( f(i) ) over all i  and n = 3-s , respectively (also see Kafatos et al. 1977). For example, in the case of codon TTA (Leu), f(1) = 1/3 (T -> C), f(2) = 0, and f(3) = 1/3 (A -> G). Thus, s = 2/3 and n = 7/3. For a DNA sequence of r codons, the total number of synonymous and nonsynonymous sites is therefore given by S = sum over j (S(i)) and n = (3r - S), respectively, where s(i) the value of s for the ith codon. When two sequences are compared, the averages of S and N for the two sequences are used.

This code implements this concept, as inspired by
https://www.biostars.org/p/60902/
and
https://github.com/a1ultima/hpcleap_dnds/blob/master/py/scripts/changes.py

inputs:
- a gff with CDS annotations
- a gff filtered for the CDS that are covered at a min read depth (get this as downfiltered from the full GFF with bedtools)
- a fasta reference genome

How does this script handle codons with "N" or "-" characters? They are ignored!
=> missing sites count as nonsyn.


"""

import argparse, gzip


# parses command line arguments
def get_commandline_arguments ():
	
	parser = argparse.ArgumentParser()
	
	parser.add_argument("--fasta", required=True, help="name/path of fasta genome reference", metavar="FILE")
	parser.add_argument("--full_gff", required=True, help="name/path of GFF annotation file, MUST be complete in terms of CDS, MUST NOT lack any CDS for any included genes", metavar="FILE")
	parser.add_argument("--filtered_gff", required=True, help="name/path of filtered GFF annotation file, MUST contain only those CDS to be included", metavar="FILE")
		
	args = parser.parse_args()

	return args



def reverse_complement (inseq):
	
	compdict = {"T":"A","C":"G","A":"T","G":"C","N":"N","-":"-"}
	
	outseq = ""
	for i in inseq[::-1]:
		outseq += compdict[i]
	
	return outseq


def read_fasta (infile):
	
	if infile.endswith(".gz"):
		F = gzip.open(infile)
	else:
		F = open(infile)
	
	fasta_dict = {}
	seq = ""
	name = "dummy"
	for line in F:
		if line.startswith(">"):
			fasta_dict[name] = seq
			name = line.lstrip(">").rstrip("\n").strip()
			print name
			seq = ""
		else:
			seq += line.strip("\n")
	
	F.close()
	
	# last record:
	fasta_dict[name] = seq
	del fasta_dict["dummy"]
	return fasta_dict


def parse_gff (infile):
	
	## NOTE: GFF coordinates are 1-based, but python indexes are 0-based
	
	## NOTE: GFF entries within each chromosome MUST be sorted by start coordinate from lowest to highest !!
	
	## GFF orientation:	+ == as-is of the fasta sequence
	##					- == reverse complement of the fasta sequence
	
	## the CDS per gene must be stitched together to interpret the codons,
	## because sometimes codons are actually split up into different exons, so that the ORF is broken when translating the split codons (frameshift)
	
	if infile.endswith(".gz"):
		F = gzip.open(infile)
	else:
		F = open(infile)
	
	gff_dict = {}
	for line in F:
		if not line.startswith("#"):
			fields = line.strip("\n").split("\t")
			if fields[2] == "CDS":
				parent = fields[8].split("Parent=")[1]
				try:
					gff_dict[parent].append(fields)
				except KeyError:
					gff_dict[parent] = [fields]
		
	return gff_dict
				

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




def interpret_gff_with_region_filter (gff_dict, fastadict, regions_to_retain, gene_whitelist):

	# purpose of gene_whitelist: there are some genes whose CDSs are overlapping in the + and - strands,
	# to avoid that both get included even though only one is explicity in the regions_to_retain, enforce the gene names with this whitelist!
	
	print "counting syn sites in CDS to be retained"
		
	out_dict = {}
	
	for k,v in gff_dict.items():
		if k in gene_whitelist:
			anyhit = None
			fourfolds = 0
			complete_CDSs = ""
			abs_coords = []
			in_gene_coordinates = []
			i = 0
			retained_CDSs = []
			for fields in v:
				chrom = fields[0]
				abs_start = int(fields[3])
				abs_stop = int(fields[4])
				# print abs_start, abs_stop
				abs_coords += range(abs_start,abs_stop+1)
				cds_length = len( range(abs_start,abs_stop+1) )
				in_gene_coordinates += range(i,cds_length+i)
				i += cds_length			
	#			ingene_stop = ingene_start + (abs_stop - abs_start)
				orientation = fields[6]
						
				try:
					complete_CDSs += fastadict[chrom][abs_start-1:abs_stop]
				except KeyError:
					complete_CDSs += ""
			
				# fo calcualtion of gene length, check if CDS is in regions to retain
				hit = site_region_lookup(abs_start, regions_to_retain[chrom])
				if hit:
					retained_CDSs.append(fields)
					anyhit = True
		
			if anyhit: # proceed only if any of the CDS in this gene is within regions_to_retain
				if len(complete_CDSs) > 0: # proceed only if the gene was actually in the fasta
					# if orientation "-", reverse the absolute coordinates list and reverse-complement the sequence
					if orientation == "-":
						abs_coords = abs_coords[::-1]
						complete_CDSs = reverse_complement (complete_CDSs)
			#		print orientation, len(v), len(in_gene_coordinates), len(abs_coords)
			#		print zip(abs_coords,in_gene_coordinates), complete_CDSs	
		
					# now estimate syn_sites
					cnt = 0
					syn_fractions_list = []
					for idx in range(0, len(complete_CDSs)-(len(complete_CDSs)%3), 3):
						codon = complete_CDSs[idx:idx+3]
						if not "N" in codon and not "-" in codon:
							syn_fractions_list += count_fraction_syn_changes (codon)
						else:
							syn_fractions_list += [0.0,0.0,0.0]
							print "ignored a codon with undetermined bases", codon
				
					# return coordinates to ascending order if necessary
					if orientation == "-":
						abs_coords = abs_coords[::-1]
						syn_fractions_list = syn_fractions_list[::-1]

					for fields in retained_CDSs:
						chrom = fields[0]
						abs_start = int(fields[3])
						abs_stop = int(fields[4])
						parent = fields[8].split("Parent=")[1]
						cdslen = abs_stop - abs_start
						syn_sites = sum( syn_fractions_list[ abs_coords.index(abs_start) : abs_coords.index(abs_stop) +1] ) 
						
#						print chrom, abs_start, abs_stop, cdslen, syn_sites, syn_sites/cdslen, parent
						
	#					print syn_fractions_list[abs_coords.index(abs_start):abs_coords.index(abs_stop)+1] 
						try:
							out_dict[chrom].append( [ abs_start, abs_stop, cdslen, syn_sites, syn_sites/cdslen, parent ] )
						except KeyError:
							out_dict[chrom] = [ [ abs_start, abs_stop, cdslen, syn_sites, syn_sites/cdslen, parent ] ]
	
	return out_dict					
	
	
def count_fraction_syn_changes (codon):
	
	# find the proportion of syn mutations among all possible mutations of the codon
	
	# translate all three possible mutations of each site of the codon, score if syn or nonsyn, count the syn ones and divide by three.
	# export syn_fraction per site
	is_AA = codon_table[codon]
	
	#print "is_codon", codon, is_AA
	
	syn_fractions = []
	
	for i in range(3):
		syn_mutations = 0
		codonlist = list(codon)
		isbase = codonlist[i]
		otherbases = [x for x in ["A","C","G","T"] if not x == isbase]
		for b in otherbases:
			x = codonlist
			x[i] = b
	#		print i, "".join(x)
			alter_AA = codon_table["".join(x)]
			if alter_AA == is_AA:
				syn_mutations += 1
		syn_fractions.append( syn_mutations/3.0 )
	
	return syn_fractions

def gff_to_regions (gff_dict):
	
	# return a SORTED list of the regions (sorted ascending by start)
	
	interimdict = {}
	for k,v in gff_dict.items():
		for fields in v:
			chrom = fields[0]
			start = int(fields[3])
			stop = int(fields[4])
			try:
				interimdict[chrom].append([start,fields])
			except KeyError:
				interimdict[chrom] = [[start,fields]]
	
	outdict = {}
	out_genelists = {}
	for k,v in interimdict.items():
		a = sorted(v, key=lambda x:x[0])
		for entry in a:
			chrom = entry[1][0]
			start = int(entry[1][3])
			stop = int(entry[1][4])
 			parent = entry[1][8].split("Parent=")[1]
 			try:
 				outdict[chrom].append((start,stop))
 				out_genelists[chrom].append(parent)
 			except KeyError:
 				outdict[chrom] = [(start,stop)]
 				out_genelists[chrom] = [parent]
	 		
 	for k,v in outdict.items():
 		x = tuple(v)
 		outdict[k] = x
 		
# 	outl = []
#	for k,y in outdict.items():
#		for v in y:
# 			outl.append(k + "\t" + str(v[0]) + "\t" + str(v[1]))
 	
# 	with open("regions_to_retain.txt", "w") as O:
# 		O.write("\n".join(outl) + "\n")
 	
 	
 	return outdict, out_genelists


def site_region_lookup(value, ranges):
	
	## https://stackoverflow.com/questions/6053974/python-efficiently-check-if-integer-is-within-many-ranges
	## "binary search"
	
	left, right = 0, len(ranges)

	while left != right - 1:
		mid = left + (right - left) // 2

		if value <= ranges[mid - 1][1]:  # Check left split max
			right = mid
		elif value >= ranges[mid][0]:    # Check right split min
			left = mid
		else:                            # We are in a gap
			return None

	if ranges[left][0] <= value <= ranges[left][1]:
		# Return the range that it fit:
#		return [ranges[left][0], ranges[left][1]]
		return True


def export_results(thedict):
	
	# within chrom, sort positions ascending, on the start coordinate
	outdict = {}
	for k,v in thedict.items():
		a = sorted(v, key=lambda x:x[0])
		outdict[k] = a
	
	outlines = [ "\t".join( ["chrom","start","stop","cdslen", "syn_sites", "syn_sites_frac", "parent_gene"] ) ]
	for k, v in outdict.items():
		for entry in v:
			outlines.append( k + "\t" + "\t".join([str(x) for x in entry]) )
	
	
#	outlines = header + [k + "\t" + "\t".join([str(x) for x in v]) for k,v in outdict.items()]
	
	with open("CDS_syn_sites.txt", "w") as O:
		O.write("\n".join(outlines) + "\n")


	
######################## MAIN

args = get_commandline_arguments ()

fastadict = read_fasta (args.fasta)

full_gff_dict = parse_gff (args.full_gff)

retainable_CDS = parse_gff (args.filtered_gff)
print "genes to retain with at least one CDS:	", len(retainable_CDS.keys())

regions_to_retain, genelists_per_chrom = gff_to_regions (retainable_CDS)

# filter directly:


chrom_dict_with_synsites_per_CDS = interpret_gff_with_region_filter(full_gff_dict, fastadict, regions_to_retain, set( retainable_CDS.keys()) )

export_results(chrom_dict_with_synsites_per_CDS)







