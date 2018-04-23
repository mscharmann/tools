#!/usr/bin/env python
"""
Unwrap fasta file so that each sequence takes up only one line.

Usage:
	%program  input_file  output_file
"""

import os, sys, argparse

# parses command line arguments
def get_commandline_arguments ():
		
	parser = argparse.ArgumentParser()
	
	parser.add_argument("--vcf", required=True, dest="vcf", type=extant_file, help="vcf genotypes", metavar="FILE")
	parser.add_argument("--ref_fasta", required=True, dest="ref_fasta", type=extant_file, help="fasta reference sequences", metavar="FILE")
	parser.add_argument("--out", required=True, help="name and path of output file", metavar="FILENAME")

	args = parser.parse_args()
	
	return args

# checks if file exists and breaks script if not
def extant_file(x):
	"""
	'Type' for argparse - checks that file exists but does not open.
	"""
	
	if not os.path.exists(x):
		print "Error: {0} does not exist".format(x)
		exit()
	x = str(x)
	return x

########

def read_wrapped_or_unwrapped_fasta (infile):
	
	outlines = []
	with open(infile, "r") as INFILE:
	
		first_id = INFILE.readline().strip("\n").strip(">")
		outlines.append(first_id)
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

def parse_vcf(vcf_file):
	
	
	with open(vcf_file, "r") as INFILE:
		for line in INFILE:
			if line.startswith("##"):
				continue
			if line.startswith("#"):
				header_line = line.lstrip("#").strip("\n").split("\t")	
				print header_line
				break
	
	samples = sorted(header_line[9:])

### freebayes vcf line:	
#	E3379_L96	43	.	ATCG	ATTG,ATCA,GTCA,GTCG	4179.3	.	AB=BLABLA	0/0:1:
	
	
### STACKS vcf line:
#	#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample_203-ampullaria-Brunei-3
#	un	2312	25	G	A	.	PASS	NS=17;AF=0.529,0.471	GT:DP:GL	1/1:32:.,.,.	

	samples_vcf_idx = {}
	for sample in samples:
		print sample
		vcf_idx = header_line.index(sample)
		samples_vcf_idx[sample] = vcf_idx
#	print popdict_vcf_idx

# now start the main business of walking through the vcf:	
	out_dict = {x:{} for x in samples}
	with open(vcf_file, "r") as INFILE:
		
		snpcnt = 0
		for line in INFILE:
			if (snpcnt / 1000.0).is_integer():
				print snpcnt
			if line.startswith("#"):
				continue
			if len(line) < 2: # empty lines or so
				continue
			fields = line.strip("\n").split("\t")
			if len(fields[3]) > 1:
				continue # exclude variants that are not SNPs; theres no way to code those things in a phylip format!
			
			snpcnt += 1
			variants = [fields[3]] + fields[4].split(",") 
			
			contig = fields[0]
			pos = fields[1]
			
			for sample in samples:
				idxes = [int(x) for x in fields[samples_vcf_idx[sample]].split(":")[0].split("/") if not x == "." ]					
				if len(idxes) > 0:
					alleles = [variants[x] for x in idxes ]
					# get IUPAC code if necessary (heterozygotes)
					outnuc = get_IUPAC_amb_code ( "".join(sorted(alleles)) ) 
				else:
					outnuc = "-"
				try:
					out_dict[sample][contig][pos] = outnuc
				except KeyError:
					out_dict[sample][contig] = {pos : outnuc}
	
		return out_dict
	

def get_IUPAC_amb_code (nucs):
	
	# nucs must be a "".joined string of the sorted list of nucleotides
	
	# the IUPAC ambiguity codes:
	# the keys in this IUPAC ambiguity code dictionary are produced by:
	# "".join(sorted(list("NUCLEOTIDES"))) -> consistent lookup possible without additional loops!
	combs = {'AC':'M', 'GT':'K', 'CG':'S', 'AT':'W', 'AG':'R', 'CT':'Y',
	'ACG':'V', 'CGT':'B', 'AGT':'D', 'ACT':'H', 'ACGT':'N'
	}

	if len(set(nucs)) > 1:
		# sometimes theres already an N, so need to except that and return N instead!
		try:		
			out_code = combs[nucs]
		except KeyError:
			out_code = "N"
	else:
		out_code = list(nucs)[0]
	return out_code


def replace_str_index(text,index=0,replacement=''):
    return '%s%s%s'%(text[:index],replacement,text[index+1:])

def genotypes_into_fasta ( refgenome_dict, samples_genotype_dict ):
	
	## we only fill the entire locus with "gap" when none of the SNPs in it was sequenced. 
	## When at least one SNP was genotxped, the locus is filled with refsequence for the other positions and 
	## with gap "-" for the SNPs that were not genotyped.
	
	gapchar = set("-")
	out_dict = {}
	for contig, positions in samples_genotype_dict.items():
		refseq = refgenome_dict[contig]
		
# 		for pos
# 		
# 		if len(refseq) < max( [ int(x) for x in positions.keys() ] ):
# 			print "FUCK!", len(refseq), [ int(x) for x in positions.keys() ]
# 			print contig, refseq
#			exit()
		
		newseq = refseq
		if len( set( positions.values() ) - gapchar ) > 0: 
			for idx, nuc in positions.items():
				if int(idx) <= len(refseq): ## need to ensure this, because apparently in some samples we have mapped reads and called SNPs on the terminal "N" that was appended to reference contigs only to allow mapping with bwa-mem; these sites are not necessarily wrong but they are undefined in samples that did not map there; so no reference contig allele exists => discard!			
#					print contig, idx, nuc
					newseq = replace_str_index(newseq, int(idx)-1, nuc)  ### we substract 1 because vcf format has 1-based offset whereas python has 0-based offset
		else:
			newseq = "-"*len( refseq )

		out_dict[contig] = newseq
		
	return out_dict


def check_for_empty_columns (concat_dict):
	
	gapchar = set("-")
	for column_index in range(0, len(concat_dict.values()[0] )): # to loop over all columns
		symbols_in_column = set( [ seq[column_index] for seq in concat_dict.values() ] )
		if len( symbols_in_column - gapchar ) == 0:
			print "dying: column ", column_index, " is empty, no output written"
			exit()
	
	print "no empty columns found, writing output"
	

				
################################## MAIN


######################

args = get_commandline_arguments ()

refgenome_file = args.ref_fasta
vcf_file = args.vcf
outfilename = args.out


ref_fasta = read_wrapped_or_unwrapped_fasta (refgenome_file)
vcf_dict = parse_vcf(vcf_file)

output_dict = {}
for sample, gts in vcf_dict.items():
	individualised_fasta_seqs = genotypes_into_fasta ( ref_fasta, gts )
	output_dict[sample] = individualised_fasta_seqs

all_contigs = vcf_dict[ vcf_dict.keys()[0] ].keys()
print "loci / contigs / alignment partitions:	", len(all_contigs)


concat_dict = {}
for i in output_dict.keys():
	models_record = []
	start_idx = 1
	end_idx = 0

	concat_seq = ""
	for l in all_contigs:
		concat_seq += output_dict[i][l]
		end_idx += len( output_dict[i][l] )
		models_record.append( "DNA," + l + "=" + str(start_idx) + "-" + str(end_idx) )
		start_idx += len( output_dict[i][l] )
	print i, len(concat_seq)
	concat_dict[i] = concat_seq

# sanity-check the alignment
check_for_empty_columns (concat_dict)

with open(args.out + ".model.txt", "w") as OUTF:
	OUTF.write( "\n".join( models_record ) + "\n")

write_phylip (concat_dict, outfilename)
print "Done!"



