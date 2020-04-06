#!/usr/local/bin/python
# Python 2.7.6
# vcf_to_genepop.py
# 3 December 2015
# Mathias Scharmann


# usage example
# python vcf_to_genepop.py --vcf bla.vcf --popfile popfile.txt

"""
- takes vcf made by freebayes and outputs a variant matrix in genepop format
- WILL NOT WORK FOR STACKS VCF DUE TO NAMING CONVENTION OF THE LOCI
- NOTE: INVARIANT sites that are called as SNPs and have not been filtered from the vcf will cause trouble downstream # so make sure the .vcf is properly filtered before.
vcftools --mac 1
- population membership is given by the popfile: one line per indiv: Ind_ID <tab> pop_ID 
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
		print "Error: {0} does not exist".format( x )
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
	parser.add_argument("--popfile", required=True, type=extant_file,
		help="name/path of population file", metavar="FILE")
	parser.add_argument("--pop_1", required=True, 
		help="name of first pop, must be same as in popfile", metavar="STRING")
	parser.add_argument("--pop_2", required=True, 
		help="name of second pop, must be same as in popfile", metavar="STRING")


# 	parser.add_argument("--format", required=True, help="the format of the vcf: stacks or freebayes")
# 	
 	args = parser.parse_args()
		
	# finish
	return args

################################## CORE

def read_pop_info (INFILE):
	
	popdict = {}
	with open(INFILE, "r") as infile:
		for line in infile:
			fields = line.strip("\n").split("\t")
			if fields[1] in popdict.keys():
				popdict[fields[1]].append(fields[0])
			else: 
				popdict[fields[1]] = [fields[0]]

	return popdict


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

	samples_vcf_idx = {}
	for sample in samples:
		print sample
		vcf_idx = header_line.index(sample)
		samples_vcf_idx[sample] = vcf_idx

	# initialise output collectors:
	genotypes_dict = {}
	for s in samples:
		genotypes_dict[s] = []
	
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
		
			for sample in samples:
				idxes = [ (int(x) + 1) for x in fields[samples_vcf_idx[sample]].split(":")[0].split("/") if not x == "." ]					
				if len(idxes) > 0:
					# we have a present genotype
					genocode = "".join( [ "0"+str(i) for i in idxes ] )						
					genotypes_dict[sample].append( genocode )
					
				else:
					# missing data
					genotypes_dict[sample].append("0000") 
	
	return genotypes_dict, loci_list


def make_genepop (genotypes_dict, loci_list, popdict, outfilename, pop_order_in_outfile):
	
	"""	
	headerline = meaningless
	lociline = all loci_IDs separated by ","
	pop
	indiv_ID,<tab>two_digits_per_allele<tab>two_digits_per_allele...
	indiv_ID,<tab>two_digits_per_allele<tab>two_digits_per_allele...
	.
	.
	.
	pop
	indiv_ID,<tab>two_digits_per_allele<tab>two_digits_per_allele...
	indiv_ID,<tab>two_digits_per_allele<tab>two_digits_per_allele...
	.
	.
	.
	"""
	headerline = "## genepop converted from .vcf for conversion to ms Hudson formats" + "	n_loci:	" + str( len(loci_list) ) + "	pops:	" + str(len( popdict.keys() )) + "	samples:	" + str( sum( [ len(x) for x in popdict.values() ] ) ) + "	pop_1 = " + pop_order_in_outfile[0] + "	pop_2 = " + pop_order_in_outfile[1]
	lociline = ",".join(loci_list)
	
	outlines = [headerline, lociline]
	for pop in pop_order_in_outfile:
		outlines.append("pop")
		for sample in popdict[pop]:
			try:
				outlines.append( sample + ",\t" + "\t".join( genotypes_dict[sample] ) )
			except KeyError:
				print "Warning: {0} on popmap but not found in vcf data".format(sample)
				continue
	
	with open(outfilename, "w") as OUTFILE:
		OUTFILE.write("\n".join( outlines ))
			
################################## MAIN


args = 	get_commandline_arguments ()
	
#print args
pop_order_in_outfile = [args.pop_1, args.pop_2]
popdict = read_pop_info(args.popfile)
for x in pop_order_in_outfile:
	if x not in popdict.keys():
		print x, " not in the popmap file"
		exit()

#print popdict
genotypes_dict, loci_list = parse_vcf(args.vcf)
print "sites:	", len(loci_list)

outfilename = args.vcf + ".genepop.txt"

make_genepop (genotypes_dict, loci_list, popdict, outfilename, pop_order_in_outfile)
	
print "Done!"
	

