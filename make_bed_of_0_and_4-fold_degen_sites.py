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
	parser.add_argument("--gtf", required=True, help="name/path of GTF annotation file, MUST be complete in terms of CDS, MUST NOT lack any CDS for any included genes", metavar="FILE")
	
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
		#	print (name)
			seq = ""
		else:
			seq += line.strip("\n")
	
	F.close()
	
	# last record:
	fasta_dict[name] = seq
	del fasta_dict["dummy"]
	return fasta_dict


def parse_gtf (infile):
	
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
				parent = fields[8].split("gene_id ")[1].split(";")[0].strip('"')
				try:
					gff_dict[parent].append(fields)
				except KeyError:
					gff_dict[parent] = [fields]
		
	print("parsed ", len(gff_dict.keys()), " genes from GTF")
	return gff_dict
				

fourfold_degen_codons = ["GCN", "CGN", "GGN", "CTN", "CCN", "TCN", "ACN", "GTN"]
# source: https://en.wikipedia.org/wiki/Codon_degeneracy
# there are 8 4-fold degenerate codons: codons that may contain any base at 1 positions (pos3) and still encode the same amino acid

fourfold_degen_codon_starts = set( [ x.rstrip("N") for x in fourfold_degen_codons ] )


# https://www.geeksforgeeks.org/dna-protein-python-3/
codon_table = {
'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',  
'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 
'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'
} 


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


def check_for_zerofold_degenerate (codon):

	# Check each position in the codon
	zero_fold_sites = []
	for pos in range(3):
		original_base = codon[pos]
		try:
			original_aa = codon_table[ codon ]
		except KeyError:
			# rare broken codons
			continue

		is_zero_fold = True

		# Test all possible substitutions
		for base in 'ACGT':
			if base == original_base:
				continue
			mutated_codon = codon[:pos] + base + codon[pos+1:]
			mutated_aa = codon_table[ mutated_codon ]
			# if at least one of the possible mutations will result in the same AA -> not 0-fold
			if mutated_aa == original_aa: 
				is_zero_fold = False
				break

		if is_zero_fold:
			zero_fold_sites.append(pos)

	return zero_fold_sites



def interpret_gtf(gff_dict, fastadict):
		
	positions_of_4xdegen = {}
	positions_of_0xdegen = {}
	
	genes_done = 0
	
	for k,v in gff_dict.items():
		genes_done += 1
		if (genes_done % 1000) == 0:
			print(genes_done)
		#print(k)
		fourfolds = 0
		complete_CDSs = ""
		abs_coords = []
		in_gene_coordinates = []
		i = 0
		retained_CDSs = []
		for fields in v:
			#print(fields)
			chrom = fields[0]
			abs_start = int(fields[3])
			abs_stop = int(fields[4])
			# print abs_start, abs_stop
			abs_coords += range(abs_start,abs_stop+1)
			cds_length = len( range(abs_start,abs_stop+1) )
			in_gene_coordinates += range(i,cds_length+i)
			i += cds_length			
	#		ingene_stop = ingene_start + (abs_stop - abs_start)
			orientation = fields[6]
						
			try:
				complete_CDSs += fastadict[chrom][abs_start-1:abs_stop]
			except KeyError:
				complete_CDSs += ""
		
		# if orientation "-", reverse the absolute coordinates list and reverse-complement the sequence
			if orientation == "-":
				abs_coords = abs_coords[::-1]
				complete_CDSs = reverse_complement (complete_CDSs)
	#		print orientation, len(v), len(in_gene_coordinates), len(abs_coords)
	#		print zip(abs_coords,in_gene_coordinates), complete_CDSs	
	
#		print("#####")
#		print(complete_CDSs)
#		print(abs_coords)
#		print(in_gene_coordinates)
		
		# now check for 4xdegen and 0xdegen
		cnt = 0
		for i in range(0, len(complete_CDSs)-(len(complete_CDSs)%3), 3):
			codon = complete_CDSs[i:i+3]
			if codon[:2] in fourfold_degen_codon_starts:
				cnt += 1
				relative_coord_4xdegen = i+2 # +1 for the 0-based index of python
				abs_coord_4xdegen = abs_coords[relative_coord_4xdegen]	
				fourfolds += 1
				try:
					positions_of_4xdegen[chrom].append(abs_coord_4xdegen)
				except KeyError:
					positions_of_4xdegen[chrom] = [abs_coord_4xdegen]
				#print("4-fold: ", chrom, abs_coord_4xdegen)
			else:
				zerofold_pos_in_codon = check_for_zerofold_degenerate (codon)
				if len( zerofold_pos_in_codon ) > 0:
					relative_coords_0xdegen = [s + i for s in zerofold_pos_in_codon]
					for s in relative_coords_0xdegen:
						abs_coord_0xdegen = abs_coords[s]
						#print("0-fold: ", chrom, abs_coord_0xdegen)
						try:
							positions_of_0xdegen[chrom].append(abs_coord_0xdegen)
						except KeyError:
							positions_of_0xdegen[chrom] = [abs_coord_0xdegen]

		
	return positions_of_4xdegen, positions_of_0xdegen
	


codon_table = {
'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',  
'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 
'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'
} 


	
######################## MAIN

args = get_commandline_arguments ()

fastadict = read_fasta (args.fasta)

full_gtf_dict = parse_gtf (args.gtf)

# filter directly:
positions_of_4xdegen, positions_of_0xdegen = interpret_gtf(full_gtf_dict, fastadict)

##

outlines = []
cnt = 0
for chrom in positions_of_4xdegen.keys():
	print ( chrom )
	for pos in positions_of_4xdegen[chrom]:
		outlines.append(chrom + "\t" + str(pos-1) + "\t" + str(pos))
		cnt += 1
	
print ("total four-fold sites:	", str(cnt))
with open("four_fold_degenerate_sites.bed", "w") as O:
	O.write("\n".join(outlines))


outlines = []
cnt = 0
for chrom in positions_of_0xdegen.keys():
	print ( chrom )
	for pos in positions_of_0xdegen[chrom]:
		outlines.append(chrom + "\t" + str(pos-1) + "\t" + str(pos))
		cnt += 1
	
print ("total zero-fold sites:	", str(cnt))
with open("zero_fold_degenerate_sites.bed", "w") as O:
	O.write("\n".join(outlines))

exit()





