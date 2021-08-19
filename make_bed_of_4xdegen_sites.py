"""
- a gff with CDS annotations
- a bed with reagions that are covered at a min read depth
- VCF with bi-allelic SNPs
- a fasta reference genome

pseudocode

before:
bedtools intersect -a -b > 
intersect GFF features with bed

- use the GFF and the FASTA to 


"""

import argparse, gzip


# parses command line arguments
def get_commandline_arguments ():
	
	parser = argparse.ArgumentParser()
	
	parser.add_argument("--fasta", required=True, help="name/path of fasta genome reference", metavar="FILE")
	parser.add_argument("--full_gff", required=True, help="name/path of GFF annotation file, MUST be complete in terms of CDS, MUST NOT lack any CDS for any included genes", metavar="FILE")
	parser.add_argument("--filtered_gff", required=True, help="name/path of filtered GFF annotation file, MUST contain only those CDS to be included", metavar="FILE")
	parser.add_argument("--vcf", required=True, help="name/path of VCF genotype file", metavar="FILE")
	
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



def interpret_gff (gff_dict, fastadict):
	
	### must construct a toy example to logic-check this function!!!
	# 1. case = single-exon gene behaves correctly
	# 2. case = multi-exon gene with exons split on lengths multiples of 3 behaves correctly
	# 3. case = multi-exon gene with exons split on lengths NOT multiples of 3: WORKS!!
		
	positions_of_4xdegen = {}
	
	genes_lengths_4xdegens = {}
	
	for k,v in gff_dict.items():
		fourfolds = 0
		remainder_bases = ""
		remainder_prev = 0
		in_gene_coordinate = -3
		for fields in v:
			chrom = fields[0]
			start = int(fields[3])
			stop = int(fields[4])
			orientation = fields[6]			
			try:
				seq = fastadict[chrom][start-1:stop]
			except KeyError:
				seq = ""
			if orientation == "-":
				seq = reverse_complement (seq)
			
			seq = remainder_bases + seq
			# handle CDS that are NOT length of multiple of 3
			remainder = int( len(seq)%3.0 )
			if remainder == 0:
				remainder_bases = ""
			else:
				remainder_bases = seq[-remainder:]
			
#			print ""
#			print k, len(seq), len(seq)%3.0, remainder_bases
#			print seq
						
			for i in range(0, len(seq)-(len(seq)%3), 3):
				in_gene_coordinate += 3
				codon = seq[i:i+3]
				if codon[:2] in fourfold_degen_codon_starts:
					fourfolds += 1
					relative_coord_4xdegen = i+2 # +1 for the 0-based index of python
					abs_coord_4xdegen = start + relative_coord_4xdegen - remainder_prev
					#print in_gene_coordinate, i, codon, "4xdegen", abs_coord_4xdegen
					try:
						positions_of_4xdegen[chrom].append(abs_coord_4xdegen)
					except KeyError:
						positions_of_4xdegen[chrom] = [abs_coord_4xdegen]
				#else:
				#	print in_gene_coordinate, i, codon
			
			remainder_prev = remainder

		# the length of the gene is defined by the start of the first exon and the end of the last exon; 
		# the coordinates are inclusive so add 1!
		genelength = int( v[-1][4] ) - int( v[0][3] ) + 1
		
		genes_lengths_4xdegens[k] = [genelength, fourfolds]		
	
	return positions_of_4xdegen, genes_lengths_4xdegens

def interpret_gff_with_region_filter_2 (gff_dict, fastadict, regions_to_retain, gene_whitelist):

	# purpose of gene_whitelist: there are some genes whose CDSs are overlapping in the + and - strands,
	# to avoid that both get included even though only one is explicity in the regions_to_retain, enforce the gene names with this whitelist!
		
	positions_of_4xdegen = {}
	
	genes_lengths_4xdegens = {}
	
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
		
				# if orientation "-", reverse the absolute coordinates list and reverse-complement the sequence
				if orientation == "-":
					abs_coords = abs_coords[::-1]
					complete_CDSs = reverse_complement (complete_CDSs)
		#		print orientation, len(v), len(in_gene_coordinates), len(abs_coords)
		#		print zip(abs_coords,in_gene_coordinates), complete_CDSs	
		
				# now check for 4xdegen
				cnt = 0
				for i in range(0, len(complete_CDSs)-(len(complete_CDSs)%3), 3):
					codon = complete_CDSs[i:i+3]
					if codon[:2] in fourfold_degen_codon_starts:
						cnt += 1
						relative_coord_4xdegen = i+2 # +1 for the 0-based index of python
						abs_coord_4xdegen = abs_coords[relative_coord_4xdegen]	
						# find out if this 4xdegen is inside of the regions to retain
						hit = site_region_lookup(abs_coord_4xdegen, regions_to_retain[chrom])
						if hit:
							#print relative_coord_4xdegen, i, codon, "4xdegen", abs_coord_4xdegen, orientation
							anyhit = True
							fourfolds += 1
							try:
								positions_of_4xdegen[chrom].append(abs_coord_4xdegen)
							except KeyError:
								positions_of_4xdegen[chrom] = [abs_coord_4xdegen]
			#	print k, "4xd sites found:	", cnt	
			
		#		if anyhit:
		#			print k, anyhit, fourfolds, orientation
				# the length of the gene is defined by the start of the first RETAINED CDS and the end of the last RETAINED CDS; 
				# the coordinates are inclusive so add 1!
				try:
					genelength = int( retained_CDSs[-1][4] ) - int( retained_CDSs[0][3] ) + 1
				except IndexError:
					genelength = 0
		
	#		print k, fourfolds, retained_CDSs, v
			if fourfolds >= 20: # genes that have less than 20 fourfolds in the regions_to_retain must be dropped.
				genes_lengths_4xdegens[k] = [genelength, fourfolds]	
	#			print chrom, k, genelength, fourfolds #, retained_CDSs	
			else:
				if anyhit:
					print "gene found in retained_regions but has no fourfold degnerate sites:	", k

		
	#	print positions_of_4xdegen
	
	return positions_of_4xdegen, genes_lengths_4xdegens
		

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

def filter_sites (sites_dict, regions_dict):
	
	out_dict = {}
	for chrom,sites in sites_dict.items():
		print "filtering 4xdegen sites to retain", chrom
		out_dict[chrom] = []
		for s in sites:
			hit = site_region_lookup(int(s), regions_dict[chrom])
			if hit:
				out_dict[chrom].append(s)

	return out_dict
		

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


def site_region_lookup_idx(value, ranges):
	
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
		# return the index of the range that the value fit in:
		return ranges.index(ranges[left])


def export_genes_lengths_4xdegens(thedict):
	
	outlines = [k + "\t" + str(v[0]) + "\t" + str(v[1]) for k,v in thedict.items()]
	
	with open("genes_lengths_4xdegens.txt", "w") as O:
		O.write("contig_name	recombination_contig_length	mutation_contig_length" + "\n")
		O.write("\n".join(outlines) + "\n")


def filter_and_modify_vcf (infile, whitelist_of_sites, regions_to_retain, genelists_per_chrom):
	
	# check if sites are among fourfolds
	# then find the gene wo which the site belongs
			
	##
	whitesets = {}
	for k, v in whitelist_of_sites.items():
		whitesets[k] = set(v)
	
	if infile.endswith(".gz"):
		INF = gzip.open(infile)
	else:
		INF = open(infile)
	
	outlines = []
	for line in INF:
		if line.startswith("#"):
			outlines.append(line.strip("\n"))
		else:
			fields = line.strip("\n").split("\t")
			chrom = fields[0]
			site = int(fields[1])
			try:
				if site in whitesets[chrom]:
					idx = site_region_lookup_idx(site,regions_to_retain[chrom])
					geneid = genelists_per_chrom[chrom][idx]
					outlines.append("\t".join( [geneid] + fields[1:] ))
					print chrom, site, geneid
			except KeyError:
				None	
	
	INF.close()
	
	with open(infile + ".4xdegen_and_chrom_as_gene.vcf", "w") as O:
		O.write("\n".join(outlines)+"\n")
	
	
######################## MAIN

args = get_commandline_arguments ()

fastadict = read_fasta (args.fasta)

full_gff_dict = parse_gff (args.full_gff)

retainable_CDS = parse_gff (args.filtered_gff)
print "genes to retain with at least one CDS:	", len(retainable_CDS.keys())

regions_to_retain, genelists_per_chrom = gff_to_regions (retainable_CDS)

# filter directly:
positions_of_4xdegen, genes_lengths_4xdegens = interpret_gff_with_region_filter_2(full_gff_dict, fastadict, regions_to_retain, set( retainable_CDS.keys()) )


outlines = []
cnt = 0
for chrom in positions_of_4xdegen.keys():
	print chrom
	for pos in positions_of_4xdegen[chrom]:
		outlines.append(chrom + "\t" + str(pos-1) + "\t" + str(pos))
		cnt += 1
	
print "total four-fold sites:	", str(cnt)
with open("four_fold_degenerate_sites.bed", "w") as O:
	O.write("\n".join(outlines))


exit()





