#!/usr/local/bin/python
# Python 2.7.6
# vcf_subset_by_chromosome.py
# 13 April 2015
# Mathias Scharmann


# usage example
# python vcf_subset_by_chromosome.py --vcf raw.01.vcf.tmp6.recode.vcf --exclude_chrom fufuchroms.txt

"""
takes vcf made by freebayes and outputs a variant matrix in phylip format for the loci subset specified in --loci

"""


import argparse
import os

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
	parser.add_argument("--exclude_chrom", required=True, type=extant_file,
		help="name/path of the file with names of chromosomes to exclude", metavar="FILE")
	
	args = parser.parse_args()
	
	# finish
	return args

################################## CORE


def parse_vcf(vcf_file, excl_chroms):
	
	try:
		os.remove(vcf_file+"subset.vcf")
	except OSError:
		None
	
	# get total lines:
	num_lines = sum(1 for line in open(vcf_file))

	# now start the main business of walking through the vcf:	
	out_lines = []
	
	with open(vcf_file, "r") as INFILE:
		
		with open(vcf_file+"subset.vcf", "a") as OUTFILE:
		
			linecnt = 0
			tlinecnt = 0
			retained_cnt = 0
			for line in INFILE:
				linecnt += 1
				tlinecnt += 1
#				print tlinecnt
				if line.startswith("#"):
					out_lines.append(line.strip("\n"))
				else:
					# throw away unwanted chromosomes:
					fields = line.strip("\n").split("\t")
					if fields[0] not in excl_chroms:
						
						# throw away also all non-biallelic SNPs:
						if "," not in fields[4]:
							out_lines.append(line.strip("\n"))
							retained_cnt += 1
					
				if linecnt == 100000:
					# write to file and flush memory
					OUTFILE.write("\n".join(out_lines)+"\n")
					linecnt = 0
					out_lines = []
				elif tlinecnt == num_lines:
					OUTFILE.write("\n".join(out_lines))
	print "retained\t" + str(retained_cnt)		
				
################################## MAIN


	
args = 	get_commandline_arguments ()
	
#print args

excl_chrom = set()
with open(args.exclude_chrom, "r") as INFILE:
	for line in INFILE:
		excl_chrom.add(line.strip("\n"))

parse_vcf(args.vcf, excl_chrom)

	
print "Done!"
	

