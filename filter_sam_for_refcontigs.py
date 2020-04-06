#!/usr/local/bin/python
# Python 2.7.6
# vcf_to_phylip.py
# 31 March 2015
# Mathias Scharmann


# usage example
# samtools view sample_TGCAT-RG.bam | python filter_sam_for_refcontigs.py --contiglist xxfu 

"""
reads STDIN from samtools view and returns only lines containing the reference contigs of names species fied in the contiglist

"""


import argparse
import os
import sys


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
	
	parser.add_argument("--contiglist", required=True, type=extant_file,
		help="name/path of contig list input file", metavar="FILE")
		
	args = parser.parse_args()
	
	# finish
	return args

################################## CORE


def read_contigfile(contiglist):
	
	contigs = set()
	with open(contiglist, "r") as INFILE:
		for line in INFILE:
			contigs.add(line.strip("\n"))
	
	return contigs


				
################################## MAIN


	
args = 	get_commandline_arguments ()

retain_contigs = read_contigfile(args.contiglist)

for line in sys.stdin:
	if line.split("\t")[2] in retain_contigs:
		sys.stdout.write(line)
	elif line.startswith("@"):
		if line.split("\t")[1].split(":")[1] in retain_contigs:
			sys.stdout.write(line)
		elif line.startswith("@RG") or line.startswith("@PG"): # retain read group in header!!
			sys.stdout.write(line)
			
#print args



	
	

