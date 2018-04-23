#!/usr/local/bin/python
# Python 2.7.6
# vcf_to_dadi.py
# 3 March 2014
# Mathias Scharmann


# usage example
# python vcf_to_dadi.py --vcf bla.vcf --popmap my_popmap.txt

"""
# example dadi script using input file prodced by this conversion:

import dadi
import pylab
dd = dadi.Misc.make_data_dict("bla.vcf.dadi_SNP_format.txt")

# for 2D SFS, find projections that maximnise fs.S() = number of segregating sites!

x = 40
res = []
for i in range(10, x):
	for j in range(10, x):
		fs = dadi.Spectrum.from_data_dict(dd , pop_ids =[ 'raff', 'amp'], projections =[i , j], polarized = False )
		s = fs.S()
		res.append([i, j, s])

		
max_S = max([x[2]for x in res])

for x in res:
	if x[2]	== max_S:
		print x

fs = dadi.Spectrum.from_data_dict(dd , pop_ids =[ 'pop1', 'pop2'], projections =[23 , 14], polarized = False )



dadi.Plotting.plot_single_2d_sfs(fs)
pylab.show()

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
	
	parser.add_argument("--popmap", required=True, type=extant_file,
		help="name and path of the population map file; format: one line per sample, sample name tabsep population name", metavar="FILE")
	parser.add_argument("--vcf", required=True, type=extant_file,
		help="name/path of vcf input file", metavar="FILE")

	
	args = parser.parse_args()
	
	linebreak_check(args.popmap)
	
	# finish
	return args

################################## CORE

def read_popmap (popmap):
	
	popdict = {}
	with open(popmap, "r") as INFILE:	
		for line in INFILE:
			sample = line.strip("\n").split("\t")[0] 
			pop = line.strip("\n").split("\t")[1] 
			try:
				popdict[pop].append(sample)
			except KeyError:
				popdict[pop] = [sample]
	return popdict

def parse_vcf(vcf_file, popdict):
	
	
	with open(vcf_file, "r") as INFILE:
		for line in INFILE:
			if line.startswith("##"):
				continue
			if line.startswith("#"):
				header_line = line.lstrip("#").strip("\n").split("\t")	
				print header_line
	
	popdict_vcf_idx = {}
	for pop in popdict.keys():
		for sample in  popdict[pop]:
			print sample
			vcf_idx = header_line.index(sample)
			try:
				popdict_vcf_idx[pop].append(vcf_idx)
			except KeyError:
				popdict_vcf_idx[pop] = [vcf_idx]
#	print popdict_vcf_idx
	
	pops = sorted(popdict.keys())
	
	with open(vcf_file, "r") as INFILE:
		with open(vcf_file+".dadi_SNP_format.txt", "w") as OUTFILE:
			
			out = ["ingroup_ref", "outgroup_ref", "Allele1"]
			for pop in pops:
				out.append(pop)
			out.append("Allele2")
			for pop in pops:
				out.append(pop)
			out.append("RADtag")
			out.append("pos")
			
			OUTFILE.write("\t".join([str(x) for x in out])+"\n")
			
			for line in INFILE:
				if line.startswith("#"):
					continue
				fields = line.strip("\n").split("\t")
				ingroup_ref = "---"
				outgroup_ref = "---"
				Allele1 = fields[3]
				Allele2 = fields[4]
				RADtag = fields[2]
				pos = fields[1] 
			
				out = [ingroup_ref, outgroup_ref, Allele1]
			
				for pop in pops:
					flds = [ fields[x] for x in popdict_vcf_idx[pop] ]
					tmp = "".join([ x.split(":")[0] for x in flds ])
					ref = tmp.count("1")
					out.append(ref)
			
				out.append(Allele2)
				for pop in pops:
					flds = [ fields[x] for x in popdict_vcf_idx[pop] ]
					tmp = "".join([ x.split(":")[0] for x in flds ])
					alt = tmp.count("0")
					out.append(alt)
			
				out.append(RADtag)
				out.append(pos)
			
				OUTFILE.write("\t".join([str(x) for x in out])+"\n")
				
################################## MAIN


	
args = 	get_commandline_arguments ()
	
#print args
	
popdict = read_popmap (args.popmap)
#print popdict

parse_vcf(args.vcf, popdict)

	
print "Done!"
	

