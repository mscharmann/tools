#!/usr/local/bin/python
# Python 2.7
# bam_to_read_depth_mock_vcf.py
# 7 Dec 2015
# Mathias Scharmann
#
# requires samtools
#
# Command line outline
# - takes no arguments

# Inputs
# bam files, should have been sorted and filtered for bwa q1 (no ambigous mappings, no unmapped reads)
	
# Outputs
# - indiv.depths.txt is the output of samtools depths (big, delete after)
# - all_mapped_sites_depths.raw.vcf
# the mock vcf can be filtered in vcftools for missingness just like a real vcf that contains useful variants
# the point is to find out the total number of observed sites, including those RADtags that had sufficient coverage but did not contain any SNPs
# Why would you want to know this? To evade ascertainment bias == overestimation of per site nucleotide diversity through the exclusion of invariant but observed sites.


import os

# get list of bam files from the current folder, each file needs to have the suffix .bam
bamfiles = [ x for x in os.listdir("./") if x.endswith(".bam")  ]

for i in bamfiles:
	print i

print "executing samtools depth"

# get the depths by calling samtools depth (may take a while!)
samtools_cmd = "samtools depth {0} > indiv.depths.txt".format( " ".join(bamfiles) )

os.system( samtools_cmd ) 

print "building mock vcf"

## build the mock vcf:
outlines = [ 
"##fileformat=VCFv4.1",
"## this is a mock-vcf just to count the individual depths on all sites",
"\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]) + "\t" + "\t".join(bamfiles) ]

cnt = 0
with open("indiv.depths.txt", "r") as INFILE:
	for line in INFILE:
		cnt += 1
		fields = line.strip("\n").split("\t")
		indiv_entries = "0/0:"+"\t0/0:".join(fields[2:])
		outline = "\t".join( [ fields[0], fields[1], ".", ".", ".", ".", ".", ".", "GT:DP", indiv_entries ])
		outlines.append(outline)

print "there are {0} mapped sites in the bam files".format(str(cnt))

with open("all_mapped_sites_depths.raw.vcf", "w") as OUTFILE:
	OUTFILE.write("\n".join(outlines))

print "Done!"

