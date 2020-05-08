#!/usr/local/bin/python
# Python 2.7.6
# vcf_select_single_random_SNP_per_RADtag.py
# 28 Januar 2016
# Mathias Scharmann


# usage example
# python vcf_select_single_random_SNP_per_RADtag.py --vcf bla.vcf

"""
- takes vcf made by freebayes and outputs a variant matrix in genepop format
- WILL NOT WORK FOR STACKS VCF DUE TO NAMING CONVENTION OF THE LOCI
- NOTE: INVARIANT sites that are called as SNPs and have not been filtered from the vcf will cause trouble downstream # so make sure the .vcf is properly filtered before.
vcftools --mac 1

- downsamples the SNPs per RADtag to have mostly independent SNPs for NeEstimator v2.01
- outputs them as list, to be filtered and then converted to genepop in different scripts.
"""


import argparse
import os
import random

########################## HEAD

# this function checks if file exists and exits script if not
def extant_file(x):
	"""
	'Type' for argparse - checks that file exists but does not open.
	"""
	
	if not os.path.exists(x):
		print "Error: {0} does not exist".format(x)
		exit()
	x = str(x)
	return x

# breaks script if non-UNIX linebreaks in input file popmap
def linebreak_check(x):

	if "\r" in open(x, "rb").readline():
		print "Error: classic mac (CR) or DOS (CRLF) linebreaks in {0}".format(x)
		exit()

######

# parses arguments from command line
def get_commandline_arguments ():
	
	parser = argparse.ArgumentParser()
	
	parser.add_argument("--vcf", required=True, type=extant_file,
		help="name/path of vcf input file", metavar="FILE")
 	args = parser.parse_args()
		
	# finish
	return args

################################## CORE

def parse_vcf(vcf_file):
	
	with open(vcf_file, "r") as INFILE:
		for line in INFILE:
			if line.startswith("##"):
				continue
			if line.startswith("#"):
				header_line = line.lstrip("#").strip("\n").split("\t")	
#				print header_line
				break
	
	samples = sorted(header_line[9:])

	### freebayes vcf line:	
	#	E3379_L96	43	.	ATCG	ATTG,ATCA,GTCA,GTCG	4179.3	.	AB=BLABLA	0/0:1:
	
	
	### STACKS vcf line:
	#	#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample_203-ampullaria-Brunei-3
	#	un	2312	25	G	A	.	PASS	NS=17;AF=0.529,0.471	GT:DP:GL	1/1:32:.,.,.	

	loci_list = []
	
	# now start the main business of walking through the vcf:		
	with open(vcf_file, "r") as INFILE:
		
		linecnt = 0
		snpcnt = 0
		for line in INFILE:
			linecnt += 1
#			print linecnt, snpcnt
			if line.startswith("#"):
				continue
			if len(line) < 2: # empty lines or so
				continue
			fields = line.strip("\n").split("\t")
			if len(fields[3]) > 1:
				continue # exclude variants that are not biallelic SNPs; we can't simulate anything else in ms anyway!
			
			# OK now we made sure that this entry is useful: NOTE: INVARIANT sites that are called as SNPs and have not been filtered from the vcf will pass this step!
			# so make sure the .vcf is properly filtered before.
			snpcnt += 1
			locusname = "-".join( [ fields[0], fields[1] ]  ) # this name format is essential for conversion script to ms Hudson
			loci_list.append( locusname )
			
	return loci_list

def sample_random_SNP (loci_list):
	
	ds_loci_list = []
	
	# decompose loci_list into RAD-tags:
	SNP_dict = {}
	# collect the idxes in the loci_list / genotypes list of SNPs per RADtag
	for idx, l in enumerate(loci_list):
		[RADtag, pos] = l.split("-")[:]
		try:
			SNP_dict[RADtag].append(idx)
		except KeyError:
			SNP_dict[RADtag] = [idx]

#	print SNP_dict

	s_loci_list = sorted(SNP_dict.keys())

	ds_loci_list = []
	for l in s_loci_list:
		picked_SNP_idx = random.randint(0, len(SNP_dict[l])-1 )
		SNP_global_idx = SNP_dict[l][picked_SNP_idx]
#		print SNP_global_idx
		ds_loci_list.append( loci_list[ SNP_global_idx ] )
			
#	print ds_genotypes_dict
#	print ds_loci_list
#	print len(ds_loci_list)

	
	return ds_loci_list

def write_loci_list (loci_list):
	
	outlines = []
	for l in loci_list:
		[RADtag, pos] = l.split("-")[:]
		outlines.append(RADtag + "\t" + pos)
	
	with open("single_SNP_sites_in_vcf.txt", "w") as OUTFILE:
		OUTFILE.write("\n".join(outlines))
		
		
################################## MAIN


	
args = 	get_commandline_arguments ()
	
loci_list = parse_vcf(args.vcf)
print "biallelic SNPs:	", len(loci_list)

loci_list = sample_random_SNP (loci_list)
print "after downsanpling to 1 random SNP per RADtag:	", len(loci_list)

write_loci_list (loci_list)
	
print "Done!"
	

