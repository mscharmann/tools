#!/usr/bin/env python
# vcf_to_PAUP_SVDQuartets_nexus.py

# Mathias Scharmann
# 2018-04-06
# python 2.7

## thanks to http://phylobotanist.blogspot.ch/2016/02/species-trees-from-snps-with-svd.html
## for explaining the PAUP SVDQuartets nexus format!

import argparse
import os
import gzip

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
		help="name/path of vcf input file, accepts gzipped format if file ends with .gz", metavar="FILE")
 	parser.add_argument("--popmap", required=True, help="the population map, tab-delim first column sample, second column pop/species", metavar="FILE")
 	parser.add_argument("--o", required=True, help="name the output file")
 	
 	args = parser.parse_args()
		
	# finish
	return args

################################## CORE


def parse_vcf(vcf_file):
	
	if vcf_file.endswith(".gz"):
		INFILE = gzip.open( vcf_file, "r" )
	else:
		INFILE = open( vcf_file, "r" )
	
	for line in INFILE:
		if line.startswith("##"):
			continue
		if line.startswith("#"):
			header_line = line.lstrip("#").strip("\n").split("\t")	
			print header_line
			break
	INFILE.close()
	samples = sorted(header_line[9:])
	samples_vcf_idx = {}
	for sample in samples:
		print sample
		vcf_idx = header_line.index(sample)
		samples_vcf_idx[sample] = vcf_idx
	samples_diploid_dict = {x:["",""] for x in samples}
# now start the main business of walking through the vcf:	
	if vcf_file.endswith(".gz"):
		INFILE = gzip.open( vcf_file, "r" )
	else:
		INFILE = open( vcf_file, "r" )
	out_columns = []
	linecnt = 0
	snpcnt = 0
	for line in INFILE:
		linecnt += 1
		if line.startswith("#"):
			continue
		if len(line) < 2: # empty lines or so
			continue
		fields = line.strip("\n").split("\t")
		if len(fields[3]) > 1:
			continue # exclude variants that are not SNPs; theres no way to code those things for SVDquartets!		
		snpcnt += 1
		column = []
		variants = [fields[3]] + fields[4].split(",") 
		# print progress every 1k sites  only
		if (snpcnt / 1000.0).is_integer(): 
			print linecnt, snpcnt
		for sample in samples:
			idxes = [int(x) for x in fields[samples_vcf_idx[sample]].split(":")[0].split("/") if not x == "." ]					
			if len(idxes) > 0:
				alleles = [variants[x] for x in idxes ]
				samples_diploid_dict[sample][0] += alleles[0]
				samples_diploid_dict[sample][1] += alleles[1]
			else:
				samples_diploid_dict[sample][0] += "?"
				samples_diploid_dict[sample][1] += "?"
	INFILE.close()
	return samples_diploid_dict


def read_popmap (infile):
	
	popdict = {}
	with open(infile, "r") as INF:
		for line in INF:
			if len( line ) > 2:
				fields = line.strip("\n").split("\t")
				try:
					popdict[ fields[1] ].append( fields[0] )
				except KeyError:
					popdict[ fields[1] ] = [ fields[0] ]
	print "read popmap"
	return popdict				



def construct_PAUP_nexus (samples_diploid_dict, popdict, outfilename):
	
	all_samples = []
	for k,v in popdict.items():
		all_samples += v
	pop_order = sorted(popdict.keys())	
	outlines = ["#nexus","begin data;"]
	ntax = len( all_samples )*2
	nchar = len( samples_diploid_dict[samples_diploid_dict.keys()[0]][0] )
	outlines.append( "   dimensions ntax = " + str(ntax) + " nchar = " + str(nchar) + " ;" )
	outlines.append( "   format datatype = dna gap = - ;" )
	outlines.append("")
	outlines.append("matrix")
	outlines.append("")
	for pop in pop_order:
		for tax in popdict[pop]:
			outlines.append("'" + tax + "_1" + "'   " + samples_diploid_dict[tax][0] )
			outlines.append("'" + tax + "_2" + "'   " + samples_diploid_dict[tax][1] )
#			outlines.append(tax + "_1" + "   " + "dummy1" )
#			outlines.append(tax + "_2" + "   " + "dummy2" )
	outlines.append("")
	outlines.append(";")
	outlines.append("end;")
	outlines.append("")
	outlines.append("begin sets;")
	outlines.append("   taxpartition yourspeciesset =")
	start = 1
	end = 0
	for pop in pop_order:
		cnt = len( popdict[pop] )*2
		end += cnt
		outlines.append("      '" + pop + "' : " + str(start) + "-" + str(end) + ",")
		start += cnt	
	tmp = outlines[-1].rstrip(",") + " ;"
	tmx = outlines[:len(outlines)-1]
#	print tmx
#	print tmp
	tmx.append(tmp)
	outlines = tmx
	outlines.append( "" )
	outlines.append( "end;" )
	with open(outfilename, "w") as OUTF:
		OUTF.write("\n".join(outlines) + "\n")
	print "Done!"


### MAIN

args = 	get_commandline_arguments ()

linebreak_check(args.popmap)

samples_diploid_dict = parse_vcf (args.vcf)
popdict = read_popmap (args.popmap)

construct_PAUP_nexus (samples_diploid_dict, popdict, args.o)


#

