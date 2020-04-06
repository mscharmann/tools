#!/usr/local/bin/python
# Python 2.7.6
# vcf_to_phylip.py
# 8 March 2014
# Mathias Scharmann


# usage example
# python vcf_to_phylip.py --vcf bla.vcf

"""
takes vcf made by freebayes and outputs a variant matrix in phylip format

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
# 	parser.add_argument("--format", required=True, help="the format of the vcf: stacks or freebayes")
# 	
 	args = parser.parse_args()
# 	
# 	if args.format not in ["stacks","freebayes"]:
# 		print "--format value not accepted; must be 'stacks' or 'freebayes' "
# 		exit()
		
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
				print header_line
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
#	print popdict_vcf_idx

# now start the main business of walking through the vcf:	
	out_columns = []
	with open(vcf_file, "r") as INFILE:
		
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
				continue # exclude variants that are not SNPs; theres no way to code those things in a phylip format!
			
			snpcnt += 1
			column = []
			variants = [fields[3]] + fields[4].split(",") 
			# print progress every 10k sites  only
			if (snpcnt / 10000.0).is_integer(): 
				print linecnt, snpcnt

			for sample in samples:
				idxes = [int(x) for x in fields[samples_vcf_idx[sample]].split(":")[0].split("/") if not x == "." ]					
				if len(idxes) > 0:
					alleles = [variants[x] for x in idxes ]
					# get IUPAC code if necessary (heterozygotes)
					column.append( get_IUPAC_amb_code ( "".join(sorted(alleles)) ) )
				else:
					column.append("-")
#			print column
			
			out_columns.append(column)
	
	# columns to rows/lines:
	outlines = []
	
	ntaxa = len(samples)
	len_align = snpcnt
	header = str(ntaxa) + " " + str(len_align)
	outlines.append(header)
	
	for idx, sample in enumerate(samples):
		seq = "".join( [ x[idx] for x in out_columns ] )
		out_line = sample + "    " + seq
		outlines.append(out_line)
	
	with open(vcf_file + ".phylip.txt", "w") as OUTFILE:
		OUTFILE.write("\n".join(outlines))

def get_IUPAC_amb_code (nucs):
	
	# nucs must be a "".joined string of the sorted list of nucleotides
	
	# the IUPAC ambiguity codes:
	# the keys in this IUPAC ambiguity code dictionary are produced by:
	# "".join(sorted(list("NUCLEOTIDES"))) -> consistent lookup possible without additional loops!
	combs = {'AC':'M', 'GT':'K', 'CG':'S', 'AT':'W', 'AG':'R', 'CT':'Y',
	'ACG':'V', 'CGT':'B', 'AGT':'D', 'ACT':'H', 'ACGT':'N'
	}

	if len(set(nucs)) > 1:
		# sometimes theres already an N, so need to except that and return N instead!
		try:		
			out_code = combs[nucs]
		except KeyError:
			out_code = "N"
	else:
		out_code = list(nucs)[0]
	return out_code
				
################################## MAIN


	
args = 	get_commandline_arguments ()
	
#print args


parse_vcf(args.vcf)

	
print "Done!"
	

