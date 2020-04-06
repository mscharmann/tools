#!/usr/local/bin/python
# Python 2.7.6
# genepop_to_mscalc_v2.5.py
# 28 April 2016
# Mathias Scharmann

"""
- converts a genepopfile for two populations into formats readible by mscalc; i.e. the spinput.txt and a genotypes file in the format of Hudson's ms coalescent simulator.
- names of SNPs in the genepopfile need to strictly follow the convention: locusID_SNPposition
- this script DOES catch uneven missingness within RADtags, which is a problem for ms Hudson and mscalc! samples with unusual missingness are dropped from that RADtag

- Allows calculation of summary statistics from observed datasets with exactly the same method as the simulation summary statistics.

- makes bpfile NOT using observed pi = theta values BUT: theta per locus = nbases*2.5*10e-8*100000
	and sample sizes per locus for simulations with msnsam

outputs:
	x.ms.txt		ready for input to mscalc
	x.spinput.txt	ready to use for calculation of observed data statistics with mscalc;
					for msnsam simulation stats calculated by mscalc using this file modify:
					 - pre-ultimate line to nreps, e.g. 10000
					 - ultimate line to name of file containing msnsam output, e.g. myfifo
	x.bpfile.txt	ready to use in msnsam simulation

example: 
python genepop_to_mscalc_v2.2.py -nreps 10000 -sp1 Hirsch -sp2 Bock -locilength 86 -gp batch_101.genepop

mv genepop_testfile.txt.spinput.txt spinput.txt

toolbase=/gdc_home3/schamath/ABC_GWH/ABC_tools/

cat genepop_testfile.txt.ms.txt | $toolbase/mscalc/mscalc

"""
import sys
import os
import argparse

# diese Funktion gibt command line arguments mit flag an das script weiter und checkt ob die Angaben plausibel / vollstandig sind
def get_commandline_arguments ():
	
	parser = argparse.ArgumentParser()
	
	parser.add_argument("-nreps", help="number of replicates / simulations")
	parser.add_argument("-sp1", help="name species 1, sp1 is first in the genepop file")
	parser.add_argument("-sp2", help="name species 2, sp2 is second in the genepop file")
#	parser.add_argument("-locilength", help="length of each locus")
	parser.add_argument("-gp", help="genepop infile", type=extant_file, required=True, metavar="FILE")
#	parser.add_argument("-obs_pi", help="tabsep file containing values of overall pi per RAD", type=extant_file, required=True, metavar="FILE")

	args = parser.parse_args()

	# finish
	return args

#######

# diese Funktion checkt ob ein file existiert und stoppt das script wenn ein file nicht exisitiert
def extant_file(x):
	"""
	'Type' for argparse - checks that file exists but does not open.
	"""
	
	if not os.path.exists(x):
		print "Error: {0} does not exist".format(x)
		exit()
	x = str(x)
	return x

###########

class Vividict(dict):
	def __missing__(self, key):
		value = self[key] = type(self)()
		return value
    
def de_vividict_loop (outer_dict1):
	
	out_dict = {}
	for key1, inner_dict1 in outer_dict1.items():
		if isinstance  (inner_dict1, Vividict):
			out_dict[key1] = dict(inner_dict1)
			outer_dict2 = inner_dict1
			for key2, inner_dict2 in outer_dict2.items():
				if isinstance  (inner_dict2, Vividict):
					out_dict[key1][key2] = dict(inner_dict2)
					outer_dict3 = inner_dict2
					for key3, inner_dict3 in outer_dict3.items():
						if isinstance  (inner_dict3, Vividict):
							out_dict[key1][key2][key3] = dict(inner_dict3)
	return out_dict	

#################


def read_genepop (genepopfile):
	
	print "reading genepop to dictionary"

	with open(genepopfile, "r") as INFILE:
		header = INFILE.readline() # remove first line
		loci_ordered = INFILE.readline().strip("\n").split(",")
		uberdict = {}
		popcount = 0
		for line in INFILE:
			items = line.strip("\n").split("\t")
			if "pop" in items:
				popcount += 1
				population = "pop_{0}".format(popcount)
				continue
			sample = "".join(items[:1]).strip(",")
			try:
				uberdict[population][sample] = items[1:]
			except KeyError:
				uberdict[population] = {sample:items[1:]}
# 				except KeyError:
# 					uberdict[population] = {sample:items[1:]}
		
	genotypes = {}
	for pop in uberdict.keys():
		genotypes[pop] = {}
		for sample in uberdict[pop].keys():
			print sample
			genotypes[pop][sample] = {}
			for idx in range(len(loci_ordered)): # previously used the zip() function here, but it is MUCH slower!
				locus = loci_ordered[idx]
				genotype = uberdict[pop][sample][idx]
#				print genotype
				genotypes[pop][sample][locus] = genotype
	
		print "read genepop to dict, now filtering"

	
	# if within an individual, the sites of a RAD-tag have unequal presence or absence, then mscalc will crash:
	# "error in reading dataset 0 : cannot read datasetfile prout"
	# this can happen when e.g. the RAD-tag is longer than the read-length and was not mapped for that individual in some positions!
	# so samples with this type of problem need be excluded altogether from that RAD -> replace with missing data the offending genotypes!
	# the error will however only be thrown if this concerns SNPs! fixed sites can have unequal missingness
	
 	replace_count = 0
 	for pop in genotypes.keys():
 		for sample in genotypes[pop].keys():
		
			locus_genotypes = [] # initiate
			locus_sites = [] # initiate
			previous_locus = loci_ordered[0].split("-")[0] # initiate
			for i in loci_ordered:
				this_locus = i.split("-")[0]
				
				if this_locus == previous_locus:
					locus_genotypes.append( genotypes[pop][sample][i] )
					locus_sites.append( i.split("-")[1]  ) # record sites of the locus
#					print locus_sites
 				else:
# 					print sample
#   					print previous_locus
#   					print locus_genotypes

 					# evaluate collected genotypes from the PREVIOUS locus:
 					pres_cnt = sum([ 1 for x in locus_genotypes if x != "0000"]) 
 					abs_cnt = locus_genotypes.count("0000")
#  					print pres_cnt
#  					print abs_cnt
 					
 					if pres_cnt < len(locus_sites):
 						if abs_cnt < len(locus_sites):
							# this means that not all SNPs in a RAD are together present or absent
												
							for site in locus_sites:
								genotypes[pop][sample]["{0}-{1}".format(previous_locus, site)] = "0000" 	
								replace_count += 1			
 					# reset collectors: but not to empty but to the values from the current loop! otherwise, loci/sites will be skipped 
 					# => the first site in each locus will be skipped, which creates problems for the replacing with missingness and hence counting later on, 
 					# since RADtag presence counter checks only the first site per RADtag...
 					locus_genotypes = [ genotypes[pop][sample][i] ] # reset collector
 					locus_sites = [ i.split("-")[1] ]				
 				
 				# update in each round
 				previous_locus = this_locus
 
 	print "warning: {0} genotypes for RAD-tags where sites within them had inconsistent number of missing data where replaced with missing data, because downstream cannot deal with this!".format(replace_count )
		
	return genotypes, loci_ordered

def genotype_stats(genotypes, loci_ordered):
	
	print "counting pres/abs of data"
	
	# this is for the spinput.txt
	loci_unique = {} # separate loci from their SNPs
	for locus in loci_ordered:
		uloc = locus.split("-")[0]
		snp = locus.split("-")[1]
		try:
			loci_unique[uloc].append(snp)
		except KeyError:
			loci_unique[uloc] = [snp]
	
	print "first stage complete"
	
	data_pres_count = {} # this works because all SNPs/sites in one locus have same number of available genotypes (or at least I thought so and will make it happen -- see filtering function above.)
	
#	sample_256-rafflesiana-Brunei-18
#	print genotypes["pop_2"]["sample_256-rafflesiana-Brunei-18"]["22304_L114-20"]
	
	for pop in genotypes.keys():
		data_pres_count[pop] = {}
		
		for locus, sites in loci_unique.items():
			data_pres_count[pop][locus] = 0
			
			for sample in genotypes[pop].keys():
				try:
					if genotypes[pop][sample]["{0}-{1}".format(locus, sites[0])] != "0000": # 0000 means missing data in this genepop format
#					print genotypes[pop][sample][locus]
						data_pres_count[pop][locus] += 2 # count 2 for each genotype because of 2 observed chromosomes
				except KeyError:
#					print "genotype missing"
					genotypes[pop][sample]["{0}-{1}".format(locus, sites[0])] = "0000"
	
	print "counted pres/abs"
					
	# bring data into a structure accessible for outputting ms-format:
	# it needs 5 dimensions...
	
	ms_datadict = {}
	for locus in loci_unique.keys():	# first level: locus
		ms_datadict[locus] = {}
		for pop in genotypes.keys():	# 2 level: populations
			ms_datadict[locus][pop] = {}
			for sample in genotypes[pop].keys():	# 3 level: samples
				ms_datadict[locus][pop][sample] = {}
				for snp in loci_unique[locus]:	# 4 level: the snps
					ms_datadict[locus][pop][sample][snp] = []	
					for chrom in [genotypes[pop][sample]["{0}-{1}".format(locus, snp)][:2], genotypes[pop][sample]["{0}-{1}".format(locus, snp)][2:]]:	# 5 level: the chromosomes, but omit missing genotypes "0000"
						if chrom != "00":
							ms_datadict[locus][pop][sample][snp].append(chrom)
	
	print "made ms_datadict with 5 dimensions"
							
	out_msdict = Vividict()	
	# now re-assign allele identifiers into 0 or 1:
	# first, make set() of alleles per snp per locus, take into account all samples and populations:
	for locus in loci_unique.keys():
		for snp in ms_datadict[locus].values()[0].values()[0].keys():
			allels = set()
			for pop in genotypes.keys():
				for sample in genotypes[pop].keys():
					for chrom in ms_datadict[locus][pop][sample][snp]:
					 allels.add(chrom)
			
			# now, replace these with 1 or 0, or omit SNP if invariant in all samples and pops
			if len(allels) == 2:
				allels_dict = {}
				allels_dict[list(allels)[0]] = "0"
				allels_dict[list(allels)[1]] = "1"
			
				for pop in genotypes.keys():
					for sample in genotypes[pop].keys():
						for chrom in [0,1]:
							# populate !
							if len(ms_datadict[locus][pop][sample][snp]) > 0: # last time getting rid of missing data/empty chrom lists
								out_msdict[locus][pop][sample][snp][chrom] = allels_dict[ ms_datadict[locus][pop][sample][snp][chrom] ] # looks up the allele and replaces allelename with ms-code

	return data_pres_count, loci_unique, out_msdict
	
def make_spinput(spinput_outfile, ms_outfile, locilength, gstats, loci_ordered):
	
	print "making spinput.txt"
	# extracts length of locus / RADtag from the name of the locus! has to be in the format 123_L100
	
	lines = []
	for locus in loci_ordered: # keep loci in order
		for pop in gstats.keys():
			lines.append(gstats[pop][locus]) # the counts of number of samples per population, multiplied by 2 for number of chromosomes!!
		locus_length = locus.split("L")[1]
		lines.append( locus_length )
		
	printlines = "\n".join([str(x) for x in lines])
	
	with open(spinput_outfile, "w") as OUTFILE: 
		OUTFILE.write("\n")
		OUTFILE.write(str(len( loci_ordered )) + "\n") # the number of loci
		OUTFILE.write(printlines + "\n")
		OUTFILE.write("1\n") # the nreps; 1 since this is only 1 observation
		OUTFILE.write(ms_outfile + "\n" + "\n") # name of the file with data in ms format
		OUTFILE.close()
	
def make_ms(ms_outfile, out_msdict, loci_ordered):

	print "making ms output"
	
	lines = []
	for locus in loci_ordered:  # keep loci in order
		lines.append("// observed locus {0}".format(locus))
		try: # some loci may not contain any segregating sites either because they were not polymorphic or because the sites were lost in filtering. Need to except that case!
			lines.append("segsites: {0}".format( len( out_msdict[locus].values()[0].values()[0].keys()) ) ) # sorry, this just finds number of seg sites
			lines.append("positions: {0}".format(" ".join([str(x) for x in  out_msdict[locus].values()[0].values()[0].keys() ]))) # this is probably meaningless
		except IndexError:
			lines.append("segsites: 0")
			lines.append("positions: -")
		
		# the core: 0 or 1 for alleles
		for pop in out_msdict[locus].keys():
			for sample in out_msdict[locus][pop].keys():
				for chrom in [0,1]:
					chromlinelist = []
#					print out_msdict[locus][pop][sample].keys()
					for snp in out_msdict[locus][pop][sample].keys():
						chromlinelist.append(out_msdict[locus][pop][sample][snp][chrom])
					chromline = "".join([str(x) for x in chromlinelist])				
					lines.append(chromline)	
				
		lines.append("")
		
	printlines = "\n".join([str(x) for x in lines])
	
	with open(ms_outfile, "w") as OUTFILE: 
		OUTFILE.write("this is observed data\n")
		OUTFILE.write("no random seeds available\n")
		OUTFILE.write("\n")
		OUTFILE.write(printlines + "\n")

	
def make_bpfile (sp1, sp2, locilength, gstats, bpoutfile, nreps, loci_ordered):
	
	print "making bpfile"
	nloci = len(loci_ordered)
	
	# give advice on msnam argument 2 
	msnsamarg = int(nloci) * int(nreps)
	print "msnsam argument 2 must be:\t{0}".format(msnsamarg)
	
	# make the output
	
	line1 = "#Sp1={0} Sp2={1}\n".format(sp1, sp2)
	
#	line2 = "\t".join([str(locilength)]*int(nloci)) + "\n"
	line2 = "\t".join( [ x.split("L")[1] for x in loci_ordered] ) + "\n"
	
	samples_pop1 = []
	for loc in loci_ordered:  # keep loci in order
		samples_pop1.append(str(gstats["pop_1"][loc]))
	line3 = "\t".join(samples_pop1) + "\n"
	
	samples_pop2 = []
	for loc in loci_ordered:  # keep loci in order
		samples_pop2.append(str(gstats["pop_2"][loc]))
	line4 = "\t".join(samples_pop2) + "\n"

	theta_values = [] # theta = observed RADlocus-specific pi
	for loc_len in [ int(x.split("L")[1]) for x in loci_ordered]:  # keep loci in order  # keep loci in order
		theta_values.append( str( loc_len* (2.5*10**(-8)) * 4 * 100000 ) )		
	line5 =	"\t".join(theta_values) + "\n"
	
	# intra-locus recombination rate rho is set to 0.0
	line6 = "\t".join(["0.0"]*int(nloci)) 
	
	with open(bpoutfile, "w") as OUTFILE: 
		OUTFILE.write(line1)
		OUTFILE.write(line2)
		OUTFILE.write(line3)
		OUTFILE.write(line4)
		OUTFILE.write(line5)
		OUTFILE.write(line6)
		OUTFILE.close()


############ MAIN

def main(argv=None):
	if argv is None:
		argv = sys.argv[1:]
		
	args = 	get_commandline_arguments ()
	
	genepopfile = args.gp

#	locilength = int(args.locilength)
	locilength = 1 # dummy value since loci length is taken directly from the observed data (the locus name contains it)
	
	spinput_outfile = "{0}.spinput.txt".format(genepopfile)
	ms_outfile = "{0}.ms.txt".format(genepopfile)
	bpoutfile = "{0}.bpfile.txt".format(genepopfile)
	
	(genotypes, loci_ordered) = read_genepop(genepopfile)
	# pop_x is the order of appearance in the genepopfile
 	
 	(gstats, loci_unique, out_msdict) = genotype_stats(genotypes, loci_ordered)
	# pop_x is still the order of appearance in the genepopfile
	
	loci_ordered = sorted( loci_unique.keys() )
	
	# all these functions produce outputs where the order of populations is the same as in the genepop input
	# sp1 is first in the genepop file, sp2 is second
 	make_bpfile (args.sp1, args.sp2, locilength, gstats, bpoutfile, args.nreps, loci_ordered)
 	
 	make_spinput(spinput_outfile, ms_outfile, locilength, gstats, loci_ordered)
 	
 	make_ms(ms_outfile, out_msdict, loci_ordered) 
		
	print "Done!"
	
	
if __name__ == "__main__":
	main(sys.argv[1:])