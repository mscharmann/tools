#!/usr/local/bin/python
# Python 2.7.6
# abba_baba_het.py M
# 02 Jan 2017
# Mathias Scharmann
# partly based on Dtest.py (pyRAD) by Deren Eaton, but IMHO much easier to input data and MUCH faster


# usage example
# python ../abba_baba_freq.17-01-02.py --i single_samples.derived_allele_freq.txt --tests chunk_${LSB_JOBINDEX} --o Dstats_out_${LSB_JOBINDEX}.txt

"""
new feature: will scan an eventual output file and only conduct the tests that have not yet been made in a previous run!

"""

import random
import itertools
import argparse, os
from math import erf, sqrt


# results / return values are:
# patterns_total = 1000 # count of all patterns analysed in all 4 taxa (roughly the number of sites)
# ABBA = 5 # count of ABBA sites
# BABA = 10 # count of BABA sites
# D = float((ABBA-BABA)) / (ABBA+BABA)
# pdisc = float((ABBA+BABA)) / patterns_total #  porportion od discordant patterns
# STD = 0.5 # standard deviation of bootstrap resampled D-statistics; bootstrapping on the loci used to calculate D
# Z = (abs(D/STD)) # the Z statistic; absolute value of D divided by the STD
# p_value = p value as converted from the Z statistic


# input: a dictionary with taxa as keys and strings of SNPs as values:
########
########

# diese Funktion gibt command line arguments mit flag an das script weiter und checkt ob die Angaben plausibel / vollstandig sind
def get_commandline_arguments ():
	
	parser = argparse.ArgumentParser()
	
	parser.add_argument("--i", required=True, type=extant_file,
		help="name and path of the derived allele frequency file", metavar="FILE")
	parser.add_argument("--tests", required=True, type=extant_file,
		help="""name and path of the tests to be done; one test per line, taxa tab-separated,\n
		p1	p2	p3	o""", metavar="FILE")
	parser.add_argument("--o", required=True, help="name and path of the output file", metavar="FILE")
	
	args = parser.parse_args()
	
	# finish
	return args

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
	

def filter_input_data (data_indict, testorder):
	# remove sites with missing data for the four taxa of this test, returns a list of lists! 
	# return also a count of the full-data sites
	
	## here we ascertain that polymorphism in {W,X} AND polymorphism in {Y,Z} because other sites make zero contribution to the D-stat (Patterson et al 2012)

	print "filtering sites for test taxa"
	indict = {}
	for tax in testorder:
		indict[tax] = data_indict[tax]
	
	testorder_freqs = []
	freqs_for_dxy = [] # here we take all sites that are present, irrespective of polymorphism in the set of taxa
	for column_index in range(0, len(indict.values()[0])): # to loop over all columns
#		print "\r" + str(column_index),
		freqs = [ indict[tax][column_index] for tax in testorder ]
		if not "-" in freqs: # require full data
			freqs = [float(x) for x in freqs]
			freqs_for_dxy.append ( freqs )
#			if not freqs[2] + freqs[3] == 0.0:
			if not freqs[2] + freqs[3] == 0.0 and not freqs[0] + freqs[1] == 0.0:
				testorder_freqs.append( freqs )
			
	
	return testorder_freqs, len(testorder_freqs), freqs_for_dxy
				

def count_informative_sites (testorder_freqs):
	
	"""
	return a count of the sites that were ABBA-BABA informative:
	sites that contain ABBA-BABA patterns MUST contain the derived allele in p3 AND contain BOTH allels in (p1,p2)
	order in testorder_freqs: [p1,p2,p3,o]	
	"""
	
	cnt = 0
	for column_index in range(0, len(testorder_freqs)): # to loop over all columns
		freqs = testorder_freqs[column_index]
		if freqs[2] > 0.0: # require derived allele presence in p3
			if sum( freqs[0:2] ) > 0.0 and sum( freqs[0:2] ) < 2.0: #sites that contain ABBA-BABA patterns MUST contain the derived allele in p3 AND contain BOTH allelea in (p1,p2)
				cnt += 1
			
	return cnt



def get_Dstat_light (patterns):

	"""
	This is applied during bootstrapping:
	The Num and Den counts for each site have already been made, so instead of doing it over and over again we just take the recorded count from prev run
	"""
		
	Num_sum = sum( [ x[0] for x in patterns])
	Den_sum = sum( [ x[1] for x in patterns])
	
	try:
		D = Num_sum / Den_sum
	except ZeroDivisionError:
		D = 0.0
	
	return D


def get_Dstat_Patterson (filtered_freq_data):
	
	"""
	order of freqs: [p1, p2, p3, O]
	
	these formulas follow:
	
	Patterson, N. J., Moorjani, P., Luo, Y., Mallick, S., Rohland, N., Zhan, Y., ... Reich, D. (2012). Ancient Admixture in Human History. Genetics, genetics.112.145037. doi:10.1534/genetics.112.145037
	
	=> hence we have here [y,z,w,x]

	ABBA = (den - num) / 2.0
	BABA = num + ((den - num) / 2.0)
	
	"""
	
	Num_sum = 0
	Den_sum = 0
	
	# illustrative only:
	ABBA_sum = 0
	BABA_sum = 0
		
	patterns = [] # record, for faster D-stat calc in resampling!
	for site in filtered_freq_data:
		# for D-stat
		y = site[0]
		z = site[1]
		w = site[2]
		x = site[3]
		
		num = (w-x)*(y-z)
		den = (w+x-2*w*x)*(y+z-2*y*z)
		
		# umgestellt aus Patterson et al p. 1702 top left
		ABBA = (den - num) / 2.0
		BABA = num + ((den - num) / 2.0)
		
		ABBA_sum += ABBA
		BABA_sum += BABA
				
		patterns.append( [ num, den ] )
		Num_sum += num
		Den_sum += den
				
	try:
		D = Num_sum / Den_sum
	except ZeroDivisionError:
		D = 0.0
	
	f_d = "na"
		
	return D, f_d, ABBA_sum, BABA_sum, patterns

		

def get_ABBA_BABA_full (filtered_freq_data):
	
	"""
	order of freqs: [p1, p2, p3, O]
	
	these formulas follow:
	
	Martin SH, Davey JW, Jiggins CD (2014) Evaluating the use of ABBA-BABA statistics to locate introgressed loci. Molecular Biology and Evolution, msu269.

	
	"""
	sum_diff_123O = 0 ## called also S(1,2,3,O)
	sum_sum_123O = 0
	sum_diff_1DDO = 0 ## called also S(1,D,D,O)
	
	# only for illustrative purpose; not actually used in calc:
	ABBA_sum = 0
	BABA_sum = 0
	
	patterns = [] # record, for faster D-stat calc in resampling!
	for site in filtered_freq_data:
		# for D-stat
		count_ABBA_123O = ( 1-site[0] )* site[1] * site[2] * ( 1-site[3] )
		count_BABA_123O = site[0] * (1-site[1]) * site[2] * ( 1-site[3] )
		
		diff_123O = count_ABBA_123O - count_BABA_123O
		sum_123O = count_ABBA_123O + count_BABA_123O
		
		patterns.append( [ diff_123O, sum_123O ] )
		ABBA_sum += count_ABBA_123O
		BABA_sum += count_BABA_123O
		
		sum_diff_123O += diff_123O
		sum_sum_123O += sum_123O
		
		# for f_d:
		if count_ABBA_123O + count_BABA_123O > 0.0:
			max23 = max( site[1], site[2] )	
			count_ABBA_1DDO = ( 1-site[0] )* max23 * max23 * ( 1-site[3] ) 
			count_BABA_1DDO = site[0] * (1-max23) * max23 * ( 1-site[3] )
		
			diff_1DDO = count_ABBA_1DDO - count_BABA_1DDO
		
			sum_diff_1DDO += diff_1DDO
	

	try:
		D = sum_diff_123O / sum_sum_123O
	except ZeroDivisionError:
		D = 0.0
	
	try:
		f_d = sum_diff_123O / sum_diff_1DDO # the proportion of the genome shared through introgression!! eq. 6 of Martin et al.
	except ZeroDivisionError:
		D = 0.0
	
	return D, f_d, ABBA_sum, BABA_sum, patterns


def sample_wr(population, k):         
    "used for bootstrap sampling"
    "Chooses k random elements (with replacement) from a population"
    "population is a range = a list of ints, k is the length of the range"
    
    n = len(population)
    _random, _int = random.random, int  # speed hack
    return [_int(_random() * n) for i in itertools.repeat(None, k)]

	
def read_tests_file (testsfile):
	
	tests = []
	with open(testsfile, "r") as INFILE:
		for line in INFILE:
			tests.append(line.strip("\n").split("\t"))
	return tests


def read_allele_freqs (alignfile):
	
	print "reading allele freqs"

	with open(alignfile, "r") as INFILE:
		pops = INFILE.readline().strip("\n").split("\t")
		data_indict = {pop : [] for pop in pops}
		for line in INFILE:
			fields = line.strip("\n").split("\t")
			for pop in pops:
				data_indict[pop].append( fields[ pops.index(pop) ] )

	return data_indict

def z2p(z):
    """
    from z-score return p-value
    """
    p=0.5*(1+erf(z/sqrt(2)))
    two_sided = 2*(1-p)
    return two_sided			

def write_results_to_file(outfilename, ret):
	
	with open(outfilename, "a") as OUTFILE:
		OUTFILE.write("\n" + "\t".join( [str(x) for x in ret ] ))

def append_header(outfilename):
	
	drawing = """
# unrooted trees:
#      /\              /\  
#     /  \   - or -   /\ \  
#    /\  /\          /\ \ \  
#   Y Z  W X        Y Z  W X  
#
# if polymorphism ascertained in {W,X} we expect D=0.  
# if D != 0 then (W,X) and (Y,Z) are not clades in the population tree
# if D < 0 then BABA-excess = gene flow between Y-W or Z-X
# if D > 0 then ABBA-excess = gene flow between Z-W or Y-X
#
"""	
	headeritems = ["Y", "Z", "W", "X", "informative_sites", "ABBA", "BABA", "D", "STD", "Z-score", "p_val", "dxy_YZ", "dxy_YW", "dxy_YX", "dxy_ZW", "dxy_ZX", "dxy_WX"]
	header = drawing + "\n" + "\t".join(headeritems)
	
	with open(outfilename, "r") as INFILE:
		fileconts = INFILE.read()
	
	if not fileconts.startswith("p1"): ## do not append header multiple times if already present!	
		with open(outfilename, "w") as OUTFILE:
			OUTFILE.write(header)
			OUTFILE.write(fileconts)
		
		
def check_already_done_tests (tests, outfilename):
		
	# which tests have already been written to the outfile?
	try:
		with open(outfilename, "r") as INF:
			done_tests = [ x.strip("\n").split("\t")[0:4] for x in INF ]	
	except IOError: # if output file does not yet exist
		done_tests = []
	
	missing_tests = [ x for x in tests if not x in done_tests ]
	print "remaining tests:	", len(missing_tests)
	
	return missing_tests


def make_dxy (all_loci_freqs):
	
	"""
	Nei & Li (1979): unnumbered eq., between eqs. 24 and 25 
	Nei 1987 eq 10.20 
	for a single biallelic SNP (1,2) in two pops it simplifies to: dxy = pop1_freq1 * ( 1.0 - pop2_freq1 ) + pop2_freq1 * ( 1.0 - pop1_freq1 )
	
	pops are ordered in inout as: [Y, Z, W, X]
	we want dxy YW and dxy ZX
	
	""" 
	dxy_YZ_list = []
	dxy_YW_list = []
	dxy_YX_list = []
	dxy_ZW_list = []
	dxy_ZX_list = []
	dxy_WX_list = []
	for s in all_loci_freqs:
#		if 0.0 < sum( s ) < 1.0: ## to exclude monomorphic sites!!
		dxy_YZ_list.append( s[0] * ( 1.0 - s[1] ) + s[1] * ( 1.0 - s[0] ) )
		dxy_YW_list.append( s[0] * ( 1.0 - s[2] ) + s[2] * ( 1.0 - s[0] ) )
		dxy_YX_list.append( s[0] * ( 1.0 - s[3] ) + s[3] * ( 1.0 - s[0] ) )
		dxy_ZW_list.append( s[1] * ( 1.0 - s[2] ) + s[2] * ( 1.0 - s[1] ) )		
		dxy_ZX_list.append( s[1] * ( 1.0 - s[3] ) + s[3] * ( 1.0 - s[1] ) )
		dxy_WX_list.append( s[3] * ( 1.0 - s[2] ) + s[2] * ( 1.0 - s[3] ) )

	try:
		dxy_YZ = sum(dxy_YZ_list) / len(dxy_YZ_list)
	except ZeroDivisionError:
		dxy_YZ = "na"		
	try:
		dxy_YW = sum(dxy_YW_list) / len(dxy_YW_list)
	except ZeroDivisionError:
		dxy_YW = "na"
	try:
		dxy_YX = sum(dxy_YX_list) / len(dxy_YX_list)
	except ZeroDivisionError:
		dxy_YX = "na"
	try:
		dxy_ZW = sum(dxy_ZW_list) / len(dxy_ZW_list)
	except ZeroDivisionError:
		dxy_ZW = "na"
	try:
		dxy_ZX = sum(dxy_ZX_list) / len(dxy_ZX_list)
	except ZeroDivisionError:
		dxy_ZX = "na"
	try:
		dxy_WX = sum(dxy_WX_list) / len(dxy_WX_list)
	except ZeroDivisionError:
		dxy_WX = "na"

	return [str( round(x, 3) ) for x in [dxy_YZ, dxy_YW, dxy_YX, dxy_ZW, dxy_ZX, dxy_WX]]
	

	
##### main #
############
############

# input: a dictionary with taxa as keys and strings of SNPs as values:

######

args = 	get_commandline_arguments ()

data_indict = read_allele_freqs (args.i)

# returns a list of lists for the different tests:
tests = read_tests_file (args.tests)

#tests = [["raff","hems","grac","perv"], ["raff","hems","grac","ghost"], ["perv","raff","grac","ghost"], ["perv","hems","grac","ghost"], ["raff","perv","grac","hems"]]

# removes
tests = check_already_done_tests (tests, args.o)

nboots = 10000

testcount = 0

for testorder in tests:
	
	testcount +=1 
	print "\n*-*-*-*-*-*-*\nD-test {0}/{1}".format(testcount, len(tests))
	
	# assign tesorder and Y Z W X:
	Y=testorder[0]
	Z=testorder[1]
	W=testorder[2]
	X=testorder[3]
	
	filtered_freq_data, count_total_sites, freqs_for_dxy = filter_input_data (data_indict, testorder)
	[dxy_YZ, dxy_YW, dxy_YX, dxy_ZW, dxy_ZX, dxy_WX] = make_dxy(freqs_for_dxy)

	obs_Dstat, obs_f_d, ABBA, BABA, ABBA_BABA_patterns = get_Dstat_Patterson (filtered_freq_data)
#	print obs_Dstat, obs_f_d, ABBA, BABA
	
	# now test significance of deviation of obs_Dstat from zero:
	print "bootstrapping"
	bs_Ds = []
	for i in range(nboots):
#		print "\r" + "bootstrapping\t" + str(i),
		# get a random list of the sites to choose:
		bs_sites = sample_wr( range( len(filtered_freq_data) ), len( filtered_freq_data ) ) 
		bs_data = [ ABBA_BABA_patterns[bs_site] for bs_site in bs_sites ]
		bs_D = get_Dstat_light (bs_data)
		bs_Ds.append( bs_D )
	
# 	with open("bs_Ds_for_hist.{}.txt".format( testcount ), "w" ) as OUTF:
# 		OUTF.write("\n".join( [str(x) for x in bs_Ds ] ))
		
	mean_bs_D = sum(bs_Ds) / len(bs_Ds)
	STD = sqrt( sum( [ (x - mean_bs_D )**2 for x in bs_Ds ]) / len(bs_Ds) )

	try:
		Zscore = (abs(obs_Dstat/STD))
	except ZeroDivisionError:
		Zscore = 0.0
		
	p_val = z2p(Zscore)
	
	ret = [Y, Z, W, X, count_total_sites, ABBA, BABA, obs_Dstat, STD, Zscore, p_val, dxy_YZ, dxy_YW, dxy_YX, dxy_ZW, dxy_ZX, dxy_WX]
	
	write_results_to_file(args.o, ret)
	
append_header (args.o)







