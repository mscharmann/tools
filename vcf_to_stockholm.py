#!/usr/local/bin/python
# Python 3
# vcf_to_stockholm.py
# 2023
# Mathias Scharmann


# usage example
# python vcf_to_stockholm.py --vcf bla.vcf

"""
takes vcf and outputs in STOCKHOLM format
https://en.wikipedia.org/wiki/Stockholm_format


this format is required by rapidnj
https://birc.au.dk/software/rapidnj

example:
rapidnj bla.vcf.sth -n -o t -a kim -c 12 -m 10000 -t d -b 100 -x bla.vcf.sth.rapidNJ_tree.tre

"""


import argparse, os, gzip

########################## HEAD

# this function checks if file exists and exits script if not
def extant_file(x):
	"""
	'Type' for argparse - checks that file exists but does not open.
	"""
	
	if not os.path.exists(x):
		print ("Error: {0} does not exist".format(x) )
		exit()
	x = str(x)
	return x


######

# parses arguments from command line
def get_commandline_arguments ():
	parser = argparse.ArgumentParser()
	parser.add_argument("--vcf", required=True, type=extant_file,
		help="name/path of vcf input file", metavar="FILE")
	args = parser.parse_args()
		
	# finish
	return args

################################## CORE


def parse_vcf(vcf_file):

	if vcf_file.endswith(".gz"):
		INFILE = gzip.open(vcf_file, "rt")
	else:
		INFILE = open(vcf_file)
	
	
	for line in INFILE:
		if line.startswith("##"):
			continue
		if line.startswith("#"):
			header_line = line.lstrip("#").strip("\n").split("\t")	
			# print (header_line)
			break
	INFILE.close()
	
	samples = sorted(header_line[9:])

### freebayes vcf line:	
#	E3379_L96	43	.	ATCG	ATTG,ATCA,GTCA,GTCG	4179.3	.	AB=BLABLA	0/0:1:
	
	
### STACKS vcf line:
#	#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample_203-ampullaria-Brunei-3
#	un	2312	25	G	A	.	PASS	NS=17;AF=0.529,0.471	GT:DP:GL	1/1:32:.,.,.	

	print ("samples in the input file: ")
	samples_vcf_idx = {}
	for sample in samples:
		print (sample)
		vcf_idx = header_line.index(sample)
		samples_vcf_idx[sample] = vcf_idx
#	print popdict_vcf_idx

# now start the main business of walking through the vcf:	
	if vcf_file.endswith(".gz"):
		INFILE = gzip.open(vcf_file, "rt")
	else:
		INFILE = open(vcf_file)

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
			continue # exclude variants that are not SNPs; theres no way to code those things in a phylip format!
			
		snpcnt += 1
#		column = []
		column = ""			
		variants = [fields[3]] + fields[4].split(",") 
		# print progress every 10k sites  only
		if (snpcnt / 10000.0).is_integer(): 
			print ("processed " + str(snpcnt) + " SNPs")
			
		for sample in samples:
			idxes = [int(x) for x in fields[samples_vcf_idx[sample]].split(":")[0].split("/") if not x == "." ]					
			if len(idxes) > 0:
				alleles = [variants[x] for x in idxes ]
				# get IUPAC code if necessary (heterozygotes)
				# column.append( get_IUPAC_amb_code ( "".join(sorted(alleles)) ) )
				column += get_IUPAC_amb_code ( "".join(sorted(alleles)) ) 
			else:
				#column.append("-")
				column += "-"
			
		out_columns.append(column)
	
	INFILE.close()
	
	# columns to rows/lines:
	outlines = []
	
	ntaxa = len(samples)
	len_align = snpcnt
	
	print("writing output:")
	print("sequences: " + str( ntaxa) )
	print("sites: " + str( len_align) )
	
	header = "# STOCKHOLM 1.0"
	outlines.append(header)
	
	
	totalgaps = 0
	for idx, sample in enumerate(samples):
		seq = "".join( [ x[idx] for x in out_columns ] )
		out_line = sample + "    " + seq
		outlines.append(out_line)
		totalgaps += seq.count("-")
		
	
	print("proportion gap: " + str( round( totalgaps/float(ntaxa*len_align), 4 ) )) 
	
	with open(vcf_file + ".sth", "w") as OUTFILE:
		OUTFILE.write("\n".join(outlines)+"\n")
		OUTFILE.write("//"+"\n")



def get_IUPAC_amb_code (nucs):
	
	# nucs must be a "".joined string of the sorted list of nucleotides
	
	# the IUPAC ambiguity codes:
	# the keys in this IUPAC ambiguity code dictionary are produced by:
	# "".join(sorted(list("NUCLEOTIDES"))) -> consistent lookup possible without additional loops!
	combs = {'AC':'M', 'GT':'K', 'CG':'S', 'AT':'W', 'AG':'R', 'CT':'Y',
	'ACG':'V', 'CGT':'B', 'AGT':'D', 'ACT':'H', 'ACGT':'N','AA':'A','TT':'T','CC':'C','GG':'G','NN':'N'
	}

	try:		
		out_code = combs[nucs]
	except KeyError:
		out_code = "N"
	
	return out_code
				
################################## MAIN


	
args = 	get_commandline_arguments ()
	
#print args


parse_vcf(args.vcf)

	
print ("Done!")
	

