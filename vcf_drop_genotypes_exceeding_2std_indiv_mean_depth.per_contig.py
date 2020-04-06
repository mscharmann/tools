#!/usr/local/bin/python
# Python 2.7.6
# vcf_drop_genotypes_exceeding_2std_indiv_mean_depth.py
# 14 April 2015
# Mathias Scharmann


# usage example
# python vcf_drop_genotypes_exceeding_2std_indiv_mean_depth.py --vcf raw.01.vcf.tmp6.recode.vcf --gdepth out.gdepth

"""
1) takes vcf made by freebayes + out.gdepth made by "vcftools --vcf fufu.vcf --geno-depth"
2) for each individual, calculates mean + 2*stdev over all sites at depth >0 (i.e. no coverage = NOT counted); gets set of RADtags (chromosomes) that exceed this threshold -> removal/blacklisted genotypes
3) parses the vcf and sets blacklisted genotypes to ".:."

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
	parser.add_argument("--gdepth", required=True, type=extant_file,
		help="name/path of the genotype depth file made with vcftools --vcf X.vcf --geno-depth", metavar="FILE")
	
	args = parser.parse_args()
	
	# finish
	return args

################################## CORE


def parse_vcf(vcf_file, tresholds):
	
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
		vcf_idx = header_line.index(sample)
		samples_vcf_idx[vcf_idx] = sample
	
	try:
		os.remove(vcf_file+".indiv_cov_filtered.vcf")
	except OSError:
		None
	
	# get total lines:
	num_lines = sum(1 for line in open(vcf_file))

	# now start the main business of walking through the vcf:	
	out_lines = []
	
	with open(vcf_file, "r") as INFILE:
		
		with open(vcf_file+".indiv_cov_filtered.vcf", "a") as OUTFILE:
		
			linecnt = 0
			tlinecnt = 0
			
			for line in INFILE:
				linecnt += 1
				tlinecnt += 1
#				print tlinecnt
				if line.startswith("#"):
					out_lines.append(line.strip("\n"))
				else:
					
					# drop blacklisted genotypes:
					fields = line.strip("\n").split("\t")
					
					outline = fields[:9]
					
					for i, record in enumerate(fields[9:]):
						vcf_idx = i + 9 
												
						try: # catch and except case where the record is only "." (e.g. as output by freebayes to indicate zero coverage)
							if not record.split(":")[1] == ".":
								dp = record.split(":")[1]
								if int(dp) <= thresholds[ samples_vcf_idx[vcf_idx] ]:
									outline.append(record)
								else:
									outline.append( "./.:" + record.split(":")[1] ) # retain only covergae information; drop GLs
									#print ".:" + ":".join(record.split(":")[1:])
							else:
								outline.append("./.:.") # retain missing data coded as ./.:. etc. 
							
						except IndexError:
							outline.append(record)
#					print "\t".join(outline)			
					out_lines.append( "\t".join(outline) )
					

					
				if linecnt == 100000:
					# write to file and flush memory
					OUTFILE.write("\n".join(out_lines)+"\n")
					linecnt = 0
					out_lines = []
				elif tlinecnt == num_lines:
					OUTFILE.write("\n".join(out_lines))
	

######

def get_thresholds (gdepth):
	
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
		
		radtag_previous = ""
		samples_max_recorder = {}
		for sample in samples:
			samples_max_recorder[sample] = 0
			
		for line in INFILE:
			linecnt += 1
#			print linecnt
						
			fields = line.strip("\n").split("\t")
			
			radtag_current = fields[0]

			if radtag_current == radtag_previous:
				
				for sample, idx in samples_idx.items():
					depth = fields[idx]
					samples_max_recorder[sample] = max([int(depth), samples_max_recorder[sample]] )

			else:
				# alter state of radtag recorder for next line:
				radtag_previous = radtag_current
				
				# record max cov of previous RADtag AND reset max recorder for the new RADtag
				for sample in samples:
					
					# the record:
					if samples_max_recorder[sample] > 0:
						samples_depths[sample].append( samples_max_recorder[sample] )
					
					samples_max_recorder[sample] = 0		
					
				for sample, idx in samples_idx.items():
					depth = fields[idx]
					samples_max_recorder[sample] = max([int(depth), 0] )	
					
#			print radtag_current, samples_max_recorder["sample_395-mirabilis-Kuching-9"]
			
#		print samples_depths
		
		thresholds = {}
		outlines = []
		for sample in samples:	
			
			depths = samples_depths[sample]
			mean = float(sum(depths))/float(len(depths))
			devs_from_mean = [ (x - mean)**2 for x in depths]
			var = sum(devs_from_mean)/len(devs_from_mean)
			sd = (var)**0.5
		
			thresholds[sample] = mean + 2*sd
			
			outlines.append("\t".join( [ str(x) for x in [ sample, len(depths), mean, var, sd, thresholds[sample] ] ] ) )
			
	
	with open( gdepth + ".indiv_coverage_report.txt", "w") as OUTFILE:
		OUTFILE.write("\t".join( ["sample", "n_RADtags_covered", "depth_mean", "var", "sd", "removal_threshold"] ) + "\n")
		OUTFILE.write("\n".join(outlines))
		
	return thresholds
#		for line in INFILE:

	
	
				
################################## MAIN


	
args = get_commandline_arguments ()
	
#print args
thresholds = get_thresholds (args.gdepth)
print thresholds

parse_vcf(args.vcf, thresholds)

	
print "Done!"
	

