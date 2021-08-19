#!/usr/local/bin/python
# Python 3
# vcf_diploid_to_haploid_random.py
# 21 January 2021
# Mathias Scharmann


# usage example
# python vcf_diploid_to_haploid_random.py --vcf bla.vcf --ploidies pl.txt

"""
takes vcf with diploid genotypes and outputs haploid genotypes, randomly sampling a single allele from each diploid genotype
=> pseudo-haploids to accomodate the Site Frequency Spectra of diploid selfers to a random-mating-like condition
"""


import argparse
import os
import random
import gzip

########################## HEAD

# this function checks if file exists and exits script if not
def extant_file(x):
	"""
	'Type' for argparse - checks that file exists but does not open.
	"""
	
	if not os.path.exists(x):
		print ("Error: {0} does not exist".format(x))
		exit()
	x = str(x)
	return x

# breaks script if non-UNIX linebreaks in input file popmap
def linebreak_check(x):

	if "\r" in open(x, "rb").readline():
		print ("Error: classic mac (CR) or DOS (CRLF) linebreaks in {0}".format(x))
		exit()

######

# parses arguments from command line
def get_commandline_arguments ():
	
	parser = argparse.ArgumentParser()
	
	parser.add_argument("--vcf", required=True, type=extant_file,
		help="name/path of vcf input file", metavar="FILE")
	parser.add_argument("--ploidies", required=True, type=extant_file,
		help="name/path of a tab-sep files with sample<tab>ploidy; allowed values are the integers 1 and 2", metavar="FILE")

	args = parser.parse_args()
# 	
# 	if args.format not in ["stacks","freebayes"]:
# 		print "--format value not accepted; must be 'stacks' or 'freebayes' "
# 		exit()
		
	# finish
	return args

################################## CORE

def parse_ploidy_file (INF):
	
	ploidy_dict = {}
	with open(INF, "r") as I:
		for line in I:
			if len(line) > 1:
				fields = line.strip().split()
				ploidy_dict[fields[0]] = fields[1]
	return ploidy_dict


def parse_vcf(vcf_file, ploidy_dict):
	
	
	# collect sample mappings
	if vcf_file.endswith(".gz"):
		INFILE = gzip.open(vcf_file, "rt")
	else:
		INFILE = open(vcf_file, "r")

	for line in INFILE:
		if line.startswith("##"):
			continue
		if line.startswith("#"):
			header_line = line.lstrip("#").strip("\n").split("\t")	
			print (header_line)
			INFILE.close()
			break
	
	samples = sorted(header_line[9:])

	samples_vcf_idx = {}
	for sample in samples:
		print (sample)
		vcf_idx = header_line.index(sample)
		samples_vcf_idx[vcf_idx] = sample

	
	# OK, now actually do the job
	outlines = []
	if vcf_file.endswith(".gz"):
		INFILE = gzip.open(vcf_file, "rt")
	else:
		INFILE = open(vcf_file, "r")
	
	snpcnt = 0
	for line in INFILE:
		if line.startswith("#"):
			outlines.append(line.strip("\n"))
			continue
		if len(line) < 2: # empty lines or so
			continue
		fields = line.strip("\n").split("\t")
		snpcnt += 1
#		print (snpcnt)
		column = []
		outl = fields[:9]
		gt_fields = fields[9:]
		field_idx = 8
		for gt in gt_fields:
			field_idx += 1
			sample = samples_vcf_idx[field_idx]
			geno = gt.split(":")[0]
			if "/" in geno:
				alleles = geno.split("/")
			elif "|" in geno:
				alleles = geno.split("|")
			if ploidy_dict[sample] == "1":
				outgeno = random.choice(alleles)
			elif ploidy_dict[sample] == "2":
				outgeno = geno
			outl.append(outgeno)
		outlines.append("\t".join(outl))

	
	INFILE.close()
	
	with open(vcf_file + ".haplo_random.vcf", "w") as O:
		O.write("\n".join(outlines)	+ "\n")	
			
				
################################## MAIN


	
args = 	get_commandline_arguments ()
	
#print args

ploidy_dict = parse_ploidy_file (args.ploidies)
print (ploidy_dict)

parse_vcf(args.vcf, ploidy_dict)

	
print ("Done!")
	

