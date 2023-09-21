#########################
# parentage_check.py
#########################

"""
Which of the candidate parents actually show highest allele sharing with a given progeny?

inputs:
- genotypes for candidate parents and progeny in VCF format
- list of candidate parent IDs
- list of progeny IDs 

output2:
- a matrix of allele sharing between progenies and each candidate parent, specifically:
	the proportion of sites where a progeny shares one or more alleles with the candidate parent
- the same matrix but instead with RANKS: ascending order from highest proportion to lowest.	

a match is defined as sharing at least one allele with the candidate parent at a given site.
a mismatch is defined as not sharing any alleles with the candidate parent at a given site.


Example usage:

python parentage_check.py --vcf test.vcf --possible_parents parents.txt --progeny progenies.txt


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
		print ("Error: {0} does not exist".format(x))
		exit()
	x = str(x)
	return x

# breaks script if non-UNIX linebreaks in input file popmap
def linebreak_check(x):

	if "\r" in open(x, "rb").readline():
		print ("Error: classic mac (CR) or DOS (CRLF) linebreaks in {0}".format(x) )
		exit()

######

# parses arguments from command line
def get_commandline_arguments ():
	
	parser = argparse.ArgumentParser()
	
	parser.add_argument("--vcf", required=True, type=extant_file,
		help="name/path of vcf input file", metavar="FILE")
	parser.add_argument("--possible_parents", required=True, type=extant_file,
		help="name/path of possible parents input file: textfile with sample names; one line per sample", metavar="FILE")
	parser.add_argument("--progeny", required=True, type=extant_file,
		help="name/path of (F1) progeny input file: textfile with sample names; one line per sample", metavar="FILE")
	args = parser.parse_args()
		
	# finish
	return args
	
#### CORE FUNCTIONS

def read_sample_list (infile):
	
	sample_list = []
	with open(infile, "r") as INFILE:
		for line in INFILE:
			if len(line) > 1:
				fields = line.strip("\n").split()
				sample_list.append(fields[0])
	return list(set(sample_list))
	

def count_allele_sharing_and_mismatch (vcf, progeny, candidate_parent):
	
	with open(vcf, "r") as INFILE:
		for line in INFILE:
			if line.startswith("##"):
				continue
			if line.startswith("#"):
				header_line = line.lstrip("#").strip("\n").split("\t")	
				break
	
	prog_i = header_line.index(progeny)
	par_i = header_line.index(candidate_parent)
	
#	print(progeny, prog_i)
#	print(candidate_parent, par_i)
	
	matches = 0
	mismatches = 0
	
	# now start the main business of walking through the vcf:	
	out_columns = []
	with open(vcf, "r") as INFILE:
		
		linecnt = 0
		snpcnt = 0
		for line in INFILE:
			linecnt += 1
			if line.startswith("#"):
				continue
			if len(line) < 2: # empty lines or so
				continue
			fields = line.strip("\n").split("\t")
			alleles_prog = set( fields[prog_i].split(":")[0].split("/") )
			alleles_par = set( fields[par_i].split(":")[0].split("/") )
			alleles_both = alleles_par.union(alleles_prog)
			# drop missing then count match or mismatch
			# print (alleles_prog, alleles_par, alleles_both)
			if not "." in alleles_both:
				if len( alleles_prog.intersection(alleles_par) ) > 0:
					matches += 1
				else: # progeny does not have at least one of the parents alleles...
					mismatches += 1
	return (round(matches /float(matches + mismatches), 10))
	



####### MAIN


	
args = 	get_commandline_arguments ()
	
#print args

potential_parents = read_sample_list(args.possible_parents)
progeny = read_sample_list(args.progeny)

print ("parents = ", potential_parents)
print ("progeny = ", progeny)


outlines1 = []
outlines1.append( "proportion of sites where a progeny shares one or more alleles with the candidate parent" )
outlines1.append( "\t".join([ "progeny_ID" ] + potential_parents) )
outlines2 = []
outlines2.append( "potential parents ranked by matching proportion, 1 = highest matching proportion" )
outlines2.append( "\t".join([ "progeny_ID" ] + potential_parents) )

for prog in progeny:
	candidate_parent_matches = []
	for par in potential_parents:
		proportion_matches = count_allele_sharing_and_mismatch (args.vcf, prog, par)
		candidate_parent_matches.append( proportion_matches )
	print(prog, candidate_parent_matches)
	outlines1.append(  "\t".join( [prog] + [str(x) for x in candidate_parent_matches ] ) )
	# now the ranks: https://stackoverflow.com/questions/3071415/efficient-method-to-calculate-the-rank-vector-of-a-list-in-python
	ranked = [sorted(candidate_parent_matches).index(x) for x in candidate_parent_matches]
	ranked1 = [abs(x-max(ranked))+1 for x in ranked] # make the rankes ascending
	outlines2.append(  "\t".join( [prog] + [str(x) for x in ranked1 ] ) )





with open("parentage_check.proportion_matches.txt", "w") as O:
	O.write("\n".join(outlines1) + "\n")



with open("parentage_check.ranks.txt", "w") as O:
	O.write("\n".join(outlines2) + "\n")


print ("Done!")
	

