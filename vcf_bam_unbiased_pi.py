#!/usr/bin/python

# 08 August 2015
# Mathias Scharmann
# vcf_bam_unbiased_pi.py

# python 2.7

# Command line outline


# usage examples


# Inputs

	
# Outputs	

#module load python/2.7 

import sys
import numpy
import os
import argparse
import subprocess

"""

"""


# checks for file existence:
def extant_file(x):
	"""
	'Type' for argparse - checks that file exists but does not open.
	"""
	
	if not os.path.exists(x):
		print "Error: {0} does not exist".format(x)
		exit()
	x = str(x)
	return x


# checks for non-UNIX linebreaks:
def linebreak_check(x):

	if "\r" in open(x, "rb").readline():
		print "Error: classic mac (CR) or DOS (CRLF) linebreaks in {0}".format(x)
		exit()
	
# parses command line arguments
def get_commandline_arguments ():
	
	parser = argparse.ArgumentParser()
	
	parser.add_argument("--bam_dir", required=True, help="name/path of vcf input file", metavar="DIRECTORY")

	parser.add_argument("--sample_list", required=True, type=extant_file,
		help="name and path of the sample list")
	
	parser.add_argument("--bam_suffix", required=True, help="suffix of the bam files, e.g. -RG.bam")
	parser.add_argument("--mindp", required=True, help="minimum depth to retain genotype", metavar="INT")
	parser.add_argument("--maxdp", required=True, help="maximum depth to retain genotype", metavar="INT")
	parser.add_argument("--max_missing", required=True, help="maximum proportion of missing genotypes to retain a locus", metavar="FLOAT")
	parser.add_argument("--vcf", required=True, type=extant_file, help="name/path of vcf input file", metavar="FILE")
	parser.add_argument("--target_contigs", required=True, type=extant_file, help="name/path of file with target contigs", metavar="FILE")
	
	args = parser.parse_args()
	
	linebreak_check(args.sample_list)
	
	return args

#######

def check_congruence (sexlistfile, bam_folder, bam_suffix):
	
	sexlistsamples = []
	with open(sexlistfile, "r") as INFILE:
		for line in INFILE:
			fields = line.strip("\n").split("\t")
			sexlistsamples.append(fields[0])
	sexlistsamples = set(sexlistsamples)
	
	for sample in sexlistsamples:
		extant_file(bam_folder + sample + bam_suffix)
	print "all samples in sex_list also have bamfiles in the bam_dir, good to go!"
	
		


#########
		 		 
def read_sexlist (sexlist_file):
	
	sexdict = {}
	with open(sexlist_file, "r") as INFILE:	
		for line in INFILE:
			sample = line.strip("\n").split("\t")[0] 
			sex = line.strip("\n").split("\t")[1] 
			if sex == "1":
				gender = "male"
			elif sex == "2":
				gender = "female"			
			try:
				sexdict[gender].append(sample)
			except KeyError:
				sexdict[gender] = [sample]
	return sexdict

###############

def parse_bam_idxstats(bam_dir, samples, bam_suffix):
	
	# read mapping info to memory; only once!
	mapping_data = {}
	for sample in samples:
		sample_mapped = get_depth_per_contig (bam_dir + sample + bam_suffix)
		for contig in sample_mapped.keys():
			try:
				mapping_data[contig].append(sample_mapped[contig])
			except KeyError:
				mapping_data[contig] = [ sample_mapped[contig] ]

	return mapping_data
	
def get_depth_per_contig (bamfile):
	
	bash_command = "samtools idxstats {0}".format(bamfile)
	print bash_command
	
	p = subprocess.Popen(bash_command, shell=True, stdout=subprocess.PIPE) #, stderr=subprocess.STDOUT)
	stdout = {}
	while True:
		line = p.stdout.readline()
#		print line,
		if line == '' and p.poll() != None:
			break
		elif "*" in line: # discard last line
			continue
		else:			
			fields = line.strip("\n").split("\t")
			# 
#			tag = fields[0]
#			dp = int(fields[2])
#			stdout[tag] = dp
			stdout[fields[0]] = int(fields[2])
	return stdout

def filter_genotype_depth (mindp, maxdp, mapping_data):
	
	# returns a dictionary that only contains loci that survive the genotype depth criteria in at least one sample
	
	filtered_data = {}
	for locus, depths in mapping_data.items():
		if sum(depths) > 0: # shortcut for completely empty loci -> dropped
			exceeds_mindp = [ x if x >= mindp else 0 for x in depths ]
			lessequal_maxdp = [ x if x <= maxdp else 0 for x in exceeds_mindp ]
			if sum(lessequal_maxdp) > 0:
				filtered_data[locus] = lessequal_maxdp
				
	return filtered_data
				
def filter_max_missing ( max_missing, in_data):
	
	min_genotypes_per_locus = int( max_missing * len(f_data[ f_data.keys()[0] ] ))
	print "discarding contigs mapped in less than {0} samples".format(min_genotypes_per_locus)
	
	out_data = {}
	for locus, depths in in_data.items():
		geno_count = len( [ x for x in depths if x != 0 ] )
		if geno_count >= min_genotypes_per_locus:
			out_data[locus] = depths

	return out_data
	
######

def parse_vcf(vcf_file, included_samples):
	
	with open(vcf_file, "r") as INFILE:
		for line in INFILE:
			if line.startswith("##"):
				continue
			if line.startswith("#"):
				header_line = line.lstrip("#").strip("\n").split("\t")	
	#			print header_line
				break
	
	#	vcf_samples = sorted(header_line[9:])

	### freebayes vcf line:	
	#	E3379_L96	43	.	ATCG	ATTG,ATCA,GTCA,GTCG	4179.3	.	AB=BLABLA	0/0:1:
	
	
	### STACKS vcf line:
	#	#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample_203-ampullaria-Brunei-3
	#	un	2312	25	G	A	.	PASS	NS=17;AF=0.529,0.471	GT:DP:GL	1/1:32:.,.,.	

	samples_vcf_idx = {}
	for sample in included_samples:
		vcf_idx = header_line.index(sample)
		samples_vcf_idx[sample] = vcf_idx

	# now start the main business of walking through the vcf:	
	out_columns = []
	with open(vcf_file, "r") as INFILE:
		
		linecnt = 0
		snpcnt = 0
		for line in INFILE:
			linecnt += 1
			print linecnt, snpcnt
			if line.startswith("#"):
				continue
			fields = line.strip("\n").split("\t")
			
			if len(fields[3]) > 1:
				continue # exclude variants that are not SNPs; theres no way to code those things in a phylip format!
			
			snpcnt += 1
			column = []
			variants = [fields[3]] + fields[4].split(",") 
		
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

####

def get_site_pi_vcftools (vcf_file, sample_list, mindp, maxdp, max_missing):
	
	bash_command = "vcftools --vcf {0} --keep {1} --mac 1 --minDP {2} --maxDP {3} --max-missing {4} --site-pi --stdout".format(vcf_file, sample_list, mindp, maxdp, max_missing)
	print bash_command
	
	p = subprocess.Popen(bash_command, shell=True, stdout=subprocess.PIPE) #, stderr=subprocess.STDOUT)
	stdout = {}
	while True:
		line = p.stdout.readline()
#		print line,
		if line == '' and p.poll() != None:
			break
		elif "CHROM" in line: # discard header
			continue
		else:			
			fields = line.strip("\n").split("\t")
			# 
#			chrom = fields[0]
#			pos = fields[1]
#			pi = float(fields[2])
#			stdout[tag] = dp
			try:
				stdout[fields[0]].append( float(fields[2]) )
			except KeyError:
				stdout[fields[0]] = [ float(fields[2]) ]
	return stdout

def combine_mapping_vcf_info (ff_data, vcf_data):
	
	# returns a dictionary { contig_ID : [ contig_length, SNP_count, [ pi_SNP1, ..., pi_SNPn ] ] }
	combined_data = {}
	for contig in ff_data.keys():
		contig_length = int(contig.split("_L")[1])
		try:
			SNP_count = len(vcf_data[contig])
			pi_list = vcf_data[contig]
		except KeyError:
			SNP_count = 0
			pi_list = []
		combined_data[contig] = [ contig_length, SNP_count, pi_list ]
	
	return combined_data
	
def make_contig_pi (combined_data, target_contigs):
	
	pi_dict = {}
	lost_contigs = set()
	for contig in target_contigs:
		try:
			pi = sum(combined_data[contig][2])/float(combined_data[contig][0])
		except KeyError:
			# if contig is not in combined_data at this point, it means that it failed the missingness or depth filters applied to the mapping data! -> simply drop it
			lost_contigs.add(contig)
			continue 
		pi_dict[contig] = pi
	
	print "total contigs targeted:	", len(target_contigs)
	print "contigs lost due to missingness / depth filtering:	", len(lost_contigs)
	print "total contigs retained:	", len(pi_dict.values())
	print "invariant contigs:	", len( [ x for x in pi_dict.values() if x == 0.0 ] )
	print "total sites seen:	", sum( [ x[0] for x in [combined_data[y] for y in target_contigs.difference(lost_contigs) ] ] )  
	print "segregating sites:	", sum( [ x[1] for x in [combined_data[y] for y in target_contigs.difference(lost_contigs) ] ] )
	print "sum of pi:	", sum( pi_dict.values() )
	print "mean pi:	",numpy.mean(pi_dict.values())
	print "var pi:	",numpy.var(pi_dict.values())
	print "std pi:	",numpy.std(pi_dict.values())
	
	print "mean pi if invariant contigs ignored:	", numpy.mean( [ x for x in pi_dict.values() if x != 0.0 ] )
	
	with open("unbiased_pi_report.txt", "w") as OUTFILE:
		outlines = []
		for contig, pi in pi_dict.items():
			outlines.append("\t".join([contig, str(pi)]))
		OUTFILE.write("\n".join(outlines))

######################## MAIN

args = get_commandline_arguments ()
	

check_congruence (args.sample_list, args.bam_dir, "-"+args.bam_suffix)

#sexdict = read_sexlist (args.samples_list)

samples = []
with open(args.sample_list) as INFILE:
	for line in INFILE:
		samp = line.strip("\n")
		samples.append(samp)
print samples

target_contigs = set()
with open(args.target_contigs) as INFILE:
	for line in INFILE:
		contig = line.strip("\n")
		target_contigs.add(contig)

####
vcf_data = get_site_pi_vcftools (args.vcf, args.sample_list, args.mindp, args.maxdp, args.max_missing)
print "contigs with biallelic SNPs:	", len(vcf_data.keys())

mapping_data = parse_bam_idxstats(args.bam_dir, samples, "-"+args.bam_suffix)
print "total contigs:	", len(mapping_data.keys())		

f_data = filter_genotype_depth (int(args.mindp), int(args.maxdp), mapping_data)
print "contigs passing depth filters:	", len(f_data.keys())

ff_data = filter_max_missing (float(args.max_missing), f_data)
print "contigs passing missingness filter:	", len(ff_data.keys())

# Now check in the .vcf: How many SNPs in the surviving contigs? If any, what is the MAF?
combined_data = combine_mapping_vcf_info (ff_data, vcf_data)

## now produce the mean pi over each contig:
make_contig_pi (combined_data, target_contigs)

#parse_mappings(args.bam_dir, samples, "-"+args.bam_suffix)

