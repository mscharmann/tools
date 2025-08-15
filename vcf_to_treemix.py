#!/usr/bin/python
# python 3
# 16 Dec 2016 etc
# Mathias Scharmann


# Command line outline


# usage example
# python 

# Inputs

	
# Outputs	

#module load python/2.7 
import os
import sys
import gzip


#######

# checks if file exists and breaks script if not
def extant_file(x):
	"""
	'Type' for argparse - checks that file exists but does not open.
	"""
	
	if not os.path.exists(x):
		print ("Error: {0} does not exist".format(x))
		exit()
	x = str(x)
	return x


# checks for non-UNIX linebreaks
def linebreak_check(x):

	if "\r" in open(x, "r").readline():
		print ("Error: classic mac (CR) or DOS (CRLF) linebreaks in {0}".format(x))
		exit()
	
# parses command line arguments
def get_commandline_arguments ():
	
	import os
	import argparse
	
	parser = argparse.ArgumentParser()
	
	parser.add_argument("--popmap", required=True, dest="popmapfile", type=extant_file, help="name and path of the popmap file: 1st column barcode/sample name separated by tab from second column indicating the population", metavar="FILE")
	parser.add_argument("--vcf", required=True, dest="vcf", type=extant_file, help="vcf genotypes", metavar="FILE")
	parser.add_argument("--max_missing", required=True, help="site-filter: maximum proportion of individuals with missing genotype per population to include a site, float value")
	
	args = parser.parse_args()
	print (args)
	
	linebreak_check(args.popmapfile)
		
	return args

########


def check_congruence (popmapfile, vcf):
	
	popmapsamples = []
	with open(popmapfile, "r") as INFILE:
		for line in INFILE:
			if len(line) > 1:
				fields = line.strip("\n").split("\t")
				popmapsamples.append(fields[0])
	popmapsamples = set(popmapsamples)
	
	
	if vcf.endswith(".gz"):
		INFILE = gzip.open(vcf, "rt")
	else:
		INFILE = open(vcf, "r")
	for line in INFILE:
		if line.startswith("##"):
			continue
		if line.startswith("#"):
			header_line = line.lstrip("#").strip("\n").split("\t")	
			break
			INFILE.close()
	
	vcf_samples = header_line[9:]
	
	for samp in popmapsamples:
		if samp not in vcf_samples:
			print ("vcf does not contain all popmap samples")
			exit()
		
	print ("all samples in popmap also in the vcf, good to go!")


#######

def read_popmap (popmapfile):
	
	popdict = {}
	with open(popmapfile,"r") as INFILE:
		for line in INFILE:
			if len(line) > 1:
				fields = line.strip("\n").split("\t")
				try:
					popdict[fields[1]].append( fields[0] )
				except KeyError:
					popdict[fields[1]] = [ fields[0] ]
	
	return popdict



def vcf_to_allele_counts (vcf, popdict, max_missing):
	
	presence_treshold = 1.0 - max_missing
	
	if vcf.endswith(".gz"):
		INFILE = gzip.open(vcf, "rt")
	else:
		INFILE = open(vcf, "r")
	for line in INFILE:
		if line.startswith("##"):
			continue
		if line.startswith("#"):
			header_line = line.lstrip("#").strip("\n").split("\t")
			break
			INFILE.close()
			
	samples = header_line[9:]
	
	#	E3379_L96	43	.	ATCG	ATTG,ATCA,GTCA,GTCG	4179.3	.	AB=BLABLA	0/0:1:
	
	popdict_vcf_idxes = {}
	for k, v in popdict.items():
		popdict_vcf_idxes[ k ] = [ header_line.index(x) for x in v ]
	
	pop_order = popdict.keys()
	
	outlines = [" ".join( pop_order )]
	
	# now start the main business of walking through the vcf:	
	print ("counting alleles")
	
	if vcf.endswith(".gz"):
		INFILE = gzip.open(vcf, "rt")
	else:
		INFILE = open(vcf, "r")

	cnt = 0
	headerlines = 0
	for line in INFILE:
		cnt += 1
		c = cnt / 10000.0
		if (c).is_integer():
			print (str(cnt))
			
		if line.startswith("#"):
			headerlines += 1				
			continue
		fields = line.strip("\n").split("\t")
		outline = []			
		for pop in pop_order:
			pop_idxs = popdict_vcf_idxes[pop]				
			a = [ fields[x].split(":")[0].split("/") for x in pop_idxs ]
			b = [x for genotype in a for x in genotype if x != "."]
			pop_n = float( len( pop_idxs ) )*2
			if len(b) >= presence_treshold * pop_n:
				p = str( b.count("1")) + "," + str( b.count("0") )
				outline.append( p )
			else:
				None
		if len(outline) == len(popdict_vcf_idxes): ## indicates that all populations survived the threshold!		
			# also make sure that this position is actually a SNP / not monomorphic in this set of populations (as suggested by M Matschiner for the F4 statistic)				
			count_A = sum( [ int( x.split(",")[0]) for x in outline ] )
			count_B = sum( [ int( x.split(",")[1]) for x in outline ] )			
			if count_A != 0 and count_B != 0: 				
				outlines.append( " ".join( outline ) )

	
	INFILE.close()
	
	print ("retained sites:	{0} out of {1}".format( len( outlines ) -1, cnt - headerlines ))
	print ("Done!")
	
	with open(vcf + "." + str(max_missing) + ".treemix", "w") as OUTF:
		OUTF.write( "\n".join( outlines ) + "\n"  )
	
	
	
		
######################## MAIN


args = get_commandline_arguments ()
	
check_congruence (args.popmapfile, args.vcf)
	
popdict = read_popmap (args.popmapfile) 
max_missing = float( args.max_missing )

vcf_to_allele_counts (args.vcf, popdict, max_missing)

