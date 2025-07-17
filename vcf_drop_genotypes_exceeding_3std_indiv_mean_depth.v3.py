#!/usr/local/bin/python
# Python 3
# vcf_drop_genotypes_exceeding_2std_indiv_mean_depth.v3.py
# Mathias Scharmann


# usage example
# zcat VCF.gz | python vcf_drop_genotypes_exceeding_2std_indiv_mean_depth.v3.py --gdepth out.gdepth | bgzip -c > OUT.vcf.gz

"""

"""


import argparse
import os, sys
import gzip

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
		print( "Error: classic mac (CR) or DOS (CRLF) linebreaks in {0}".format(x))
		exit()

######

# parses arguments from command line
def get_commandline_arguments ():
	
	parser = argparse.ArgumentParser()

	parser.add_argument("--gdepth", required=True, type=extant_file,
		help="name/path of the genotype depth file made with vcftools --vcf X.vcf --geno-depth", metavar="FILE")
	
	args = parser.parse_args()
	
	# finish
	return args

################################## CORE




def parse_vcf(thresholds):
	
	# this was improved for speed with chatGPT..
	
	for line in sys.stdin:
		if line.startswith("##"):
			sys.stdout.write(line )
			continue

		if line.startswith("#"):
			header_line = line.lstrip("#").strip().split("\t")
			samples = sorted(header_line[9:])
			samples_vcf_idx = {header_line.index(sample): sample for sample in samples}
			surviving_counter = {sample: 0 for sample in samples}
			sys.stdout.write(line )
			continue

		# Process variant lines
		fields = line.strip().split("\t")
		outline = fields[:9]

		for i, record in enumerate(fields[9:]):
			vcf_idx = i + 9
			# Early handling of missing data
			if record == ".":
				outline.append(".")
				continue

			parts = record.split(":")
			dp = None

			try:
				dp = int(parts[1])
			except ValueError:
				try:
					dp = int(parts[2])
				except (ValueError, IndexError):
					outline.append(record)
					continue

			sample_name = samples_vcf_idx[vcf_idx]

			if dp is not None and dp <= thresholds[sample_name]:
				outline.append(record)
				surviving_counter[sample_name] += 1
			else:
				outline.append(".:" + parts[1] if len(parts) > 1 else ".")

		sys.stdout.write("\t".join(outline) + "\n")

	return surviving_counter

######

import random

def reservoir_sampling(file_path, N):
	"""
	Randomly sample N lines from a file using reservoir sampling.

	:param file_path: Path to the input file
	:param N: Number of lines to sample
	:return: List of N randomly sampled lines
	"""
	reservoir = []
	with open(file_path, 'r') as file:
		file.readline() # drop the header
		for i, line in enumerate(file):
			if i < N:
				reservoir.append(line.strip())
			else:
				# Randomly replace elements in the reservoir with a decreasing probability
				j = random.randint(0, i)
				if j < N:
					reservoir[j] = line.strip()
	return reservoir


def get_thresholds (gdepth):
	
	# randomly sample 10k lines from the file
	sampled_lines = reservoir_sampling(gdepth, 20000)
	
	samples_depths = {}
	
	with open(gdepth, "r") as INFILE:

		header = INFILE.readline().strip("\n").split("\t")
		## CHROM	POS	sample_100-hemsleyana-Brunei-4	sample_102-rafflesiana-Brunei-4	...
		samples = sorted(header[2:])
		samples_idx = {}
		for sample in samples:
			idx = header.index(sample)
			samples_idx[sample] = idx
			
			#samples_depths[sample] = [0,0] # dictionary records [count_of_covered_sites, sum_of_reads]
			samples_depths[sample] = []
		
		linecnt = 0	
					
	for line in sampled_lines:
		linecnt += 1
#		print linecnt
						
		fields = line.strip("\n").split("\t")
		

		for sample, idx in samples_idx.items():
			depth = float(fields[idx])
			samples_depths[sample].append(depth)
		
	thresholds = {}
	for sample in samples:	
			
		nonzero_depths = [x for x in samples_depths[sample] if not x == 0.0]

		try:
			mean = float(sum(nonzero_depths))/float(len(nonzero_depths))
			devs_from_mean = [ (x - mean)**2 for x in nonzero_depths]
			var = sum(devs_from_mean)/len(devs_from_mean)
			sd = (var)**0.5
		except ZeroDivisionError:
			mean = 0.0
			devs_from_mean = [ (x - mean)**2 for x in nonzero_depths]
			var = 0.0
			sd = (var)**0.5
			
		thresholds[sample] = mean + 3*sd
							
	
	#print(thresholds)		
	return thresholds, samples_depths
#		for line in INFILE:


def write_report (samples_depths, surviving_counter, gdepth):
	
	thresholds = {}
	outlines = []
	for sample in surviving_counter.keys():	
			
		nonzero_depths = [x for x in samples_depths[sample] if not x == 0.0]
		try:
			mean = float(sum(nonzero_depths))/float(len(nonzero_depths))
			devs_from_mean = [ (x - mean)**2 for x in nonzero_depths]
			var = sum(devs_from_mean)/len(devs_from_mean)
			sd = (var)**0.5
		except ZeroDivisionError:
			mean = 0.0
			devs_from_mean = [ (x - mean)**2 for x in nonzero_depths]
			var = 0.0
			sd = (var)**0.5
	
		thresholds[sample] = mean + 3*sd
		
		outlines.append("\t".join( [ str(x) for x in [ sample, len(nonzero_depths), mean, var, sd, thresholds[sample], surviving_counter[sample] ] ] ) )

	
	with open( gdepth + ".indiv_coverage_report.txt", "w") as OUTFILE:
		OUTFILE.write("\t".join( ["sample", "nonzero_sites_sampled_for_depth_estimates", "depth_mean", "var", "sd", "removal_threshold", "sites_surviving_filtering"] ) + "\n")
		OUTFILE.write("\n".join(outlines))
	
	
				
################################## MAIN


	
args = get_commandline_arguments ()
	
#print args
thresholds, samples_depths = get_thresholds (args.gdepth)
#print thresholds

surviving_counter = parse_vcf(thresholds)
write_report (samples_depths, surviving_counter, args.gdepth)

	
#print ("Done!")
	

