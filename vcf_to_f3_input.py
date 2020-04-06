#!/usr/bin/python
# python 2.7
# 15 Dec 2016
# Mathias Scharmann


# Command line outline


# usage example
# python 

# Inputs

	
# Outputs	

#module load python/2.7 
import os
import sys


#######

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


# checks for non-UNIX linebreaks
def linebreak_check(x):

	if "\r" in open(x, "rb").readline():
		print "Error: classic mac (CR) or DOS (CRLF) linebreaks in {0}".format(x)
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
	print args
	
	linebreak_check(args.popmapfile)
		
	return args

########


def check_congruence (popmapfile, vcf):
	
	popmapsamples = []
	with open(popmapfile, "r") as INFILE:
		for line in INFILE:
			fields = line.strip("\n").split("\t")
			popmapsamples.append(fields[0])
	popmapsamples = set(popmapsamples)
	
	with open(vcf, "r") as INFILE:
		for line in INFILE:
			if line.startswith("##"):
				continue
			if line.startswith("#"):
				header_line = line.lstrip("#").strip("\n").split("\t")	
				break
	
	vcf_samples = header_line[9:]
	
	for samp in popmapsamples:
		if samp not in vcf_samples:
			print "vcf does not contain all popmap samples"
			exit()
	for samp in vcf_samples:
		if samp not in popmapsamples:
			print "popmap does not contain all vcf samples"
			exit()
		
	print "samples in popmap and vcf fully congruent, good to go!"


#######

def read_popmap (popmapfile):
	
	popdict = {}
	with open(popmapfile,"r") as INFILE:
		for line in INFILE:
			fields = line.strip("\n").split("\t")
			try:
				popdict[fields[1]].append( fields[0] )
			except KeyError:
				popdict[fields[1]] = [ fields[0] ]
	
	return popdict



def vcf_to_allele_freqs (vcf, popdict, max_missing):
	
	presence_treshold = 1.0 - max_missing
	
	with open(vcf, "r") as INFILE:
		for line in INFILE:
			if line.startswith("##"):
				continue
			if line.startswith("#"):
				header_line = line.lstrip("#").strip("\n").split("\t")
				break
	
	samples = header_line[9:]
	
	#	E3379_L96	43	.	ATCG	ATTG,ATCA,GTCA,GTCG	4179.3	.	AB=BLABLA	0/0:1:
	
	popdict_vcf_idxes = {}
	for k, v in popdict.items():
		popdict_vcf_idxes[ k ] = [ header_line.index(x) for x in v ]
	
	# now start the main business of walking through the vcf:	
	print "counting allele frequencies"
	
	allele_freq_dict = {k : [] for k in popdict.keys() }
	with open(vcf, "r") as INFILE:
		
		cnt = 0
		for line in INFILE:
			cnt += 1
			print str(cnt) + "\r",
			
			if line.startswith("#"):
				continue
			fields = line.strip("\n").split("\t")
			for pop, pop_idxs in popdict_vcf_idxes.items():
				a = [ fields[x].split(":")[0].split("/") for x in pop_idxs ]
				b = [x for genotype in a for x in genotype if x != "."]
				pop_n = float( len( pop_idxs ) )*2
				if len(b) >= presence_treshold * pop_n:
					p = str( b.count("1")/float(len(b)) )
				else:
					p = "-"
				allele_freq_dict[pop].append( p )
				
	
	# clean allele_freq_dict from sites that are not present in at least three populations:
	print "filtering sites in allele frequency table for presence in >= 3 populations"
	pop_order = allele_freq_dict.keys()
	output_lines = [pop_order]
	
	for i in range( len( allele_freq_dict[allele_freq_dict.keys()[0]] ) ):
		column = [ allele_freq_dict[pop][i] for pop in pop_order ]
		if len( [ x for x in column if x != "-" ] ) >= 3:
			output_lines.append( column ) # now its a row!

	print "retained sites:	{0} out of {1}".format( len( output_lines ) -1, len( allele_freq_dict[allele_freq_dict.keys()[0]] ) )
	print "Done!"
	
	with open(vcf + ".maxmiss_per_pop." + str(max_missing) + ".allele_freqs.txt", "w") as OUTF:
		OUTF.write( "\n".join(["\t".join(x) for x in output_lines ])  )
	
	
	
		
######################## MAIN


args = get_commandline_arguments ()
	
check_congruence (args.popmapfile, args.vcf)
	
popdict = read_popmap (args.popmapfile) 
max_missing = float( args.max_missing )

vcf_to_allele_freqs (args.vcf, popdict, max_missing)

