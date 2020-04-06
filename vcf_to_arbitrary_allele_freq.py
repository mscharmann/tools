#!/usr/bin/python
# python 2.7
# 21 Dec 2016
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
	parser.add_argument("--max_missing_ingr", required=True, help="site-filter: maximum proportion of individuals with missing genotype per population to include a site, float value")
	parser.add_argument("--out", required=True, help="name and path of output file", metavar="FILENAME")

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
		
	print "all samples in popmap are in .vcf, good to go!"


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



def vcf_to_arbitrary_allele_freqs (vcf, popdict, max_missing_ingr, output_file):
	
	presence_treshold_ingr = 1.0 - max_missing_ingr
	
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
	
	ingroup_pops = [ x for x in popdict.keys() ]
	
	pop_order = popdict.keys()
	
	output_lines = [pop_order]
	with open(output_file, "w") as OUTF:
		OUTF.write( "\n".join(["\t".join(x) for x in output_lines ]) + "\n"  )
	
	# now start the main business of walking through the vcf:	
	print "counting allele frequencies, polarisation is arbitrary (NOT necessarily minor or major.)"
	
	allele_freq_dict = {k : [] for k in pop_order }

	with open(vcf, "r") as INFILE:
		
		cnt = 0
		for line in INFILE:
			if line.startswith("#"):
				continue

			cnt += 1
			c = cnt / 100000.0
			if (c).is_integer():
				print str(cnt)
				write_ouput_chunk ( allele_freq_dict, pop_order, output_file )
				# and clean memory by replacing the dict:
				allele_freq_dict = {k : [] for k in pop_order }

			fields = line.strip("\n").split("\t")
			## first check that site occurs in the outgroup and fixed in the outgroup:
#			a = [ fields[x].split(":")[0].split("/") for x in popdict_vcf_idxes[outgr] ]
#			b = [x for genotype in a for x in genotype if x != "."]
#			pop_n = float( len( popdict_vcf_idxes[outgr] ) )*2
#			if len(b) >= presence_treshold_outgr * pop_n:
#				if len(set(b)) == 1:
					# we have a site that survives max_missing AND is fixed in the outgroup
					# now check the other populations:
#			print fields
			for pop in ingroup_pops:
				a = [ fields[x].split(":")[0].split("/") for x in popdict_vcf_idxes[pop] ]
				b = [x for genotype in a for x in genotype if x != "."]
				pop_n = float( len( popdict_vcf_idxes[pop] ) )*2
#				print pop_n
				if len(b) >= presence_treshold_ingr * pop_n:
					p = str( b.count( "1" )/float(len(b)) )   # OK we just count the allele "1", it does not matter
				else:
					p = "-"
				allele_freq_dict[pop].append( p )
				
	# final write:
	write_ouput_chunk ( allele_freq_dict, pop_order, output_file )	
	return cnt


def write_ouput_chunk ( allele_freq_dict, pop_order, output_file ):

	"""
	better to write the output from time to time rather than to clog up the memory with GB! 
	"""
	
	# clean allele_freq_dict from sites that are not present in at least four populations:
	print "filtering sites in allele frequency table for presence in >= 4 populations"
	
	output_lines = []
	for i in range( len( allele_freq_dict[allele_freq_dict.keys()[0]] ) ):
		column = [ allele_freq_dict[pop][i] for pop in pop_order ]
		if len( [ x for x in column if x != "-" ] ) >= 4:
			output_lines.append( column ) # now its a row!

	with open(output_file, "a") as OUTF:
		OUTF.write( "\n".join(["\t".join(x) for x in output_lines ]) + "\n"  )

		
######################## MAIN


args = get_commandline_arguments ()
	
check_congruence (args.popmapfile, args.vcf)
	
popdict = read_popmap (args.popmapfile) 
max_missing_ingr = float( args.max_missing_ingr )

input_sites_cnt = vcf_to_arbitrary_allele_freqs (args.vcf, popdict, max_missing_ingr, args.out)

## check length of output_file:
output_lines = sum(1 for line in open( args.out )) -1

print "retained sites:	{0} out of {1}".format( output_lines, input_sites_cnt )
print "Done!"


