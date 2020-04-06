#!/usr/bin/python
# python 2.7
# 09 Oct 2017
# Mathias Scharmann


# Command line outline


# usage example
# python 

# Inputs

	
# Outputs	

#module load python/2.7 
import os
import sys


#######

# checks if file exists and breaks script if not
def extant_file(x):
	"""
	'Type' for argparse - checks that file exists but does not open.
	"""
	
	if not os.path.exists(x):
		print "Error: {0} does not exist".format(x)
		exit()
	x = str(x)
	return x


# checks for non-UNIX linebreaks
def linebreak_check(x):

	if "\r" in open(x, "rb").readline():
		print "Error: classic mac (CR) or DOS (CRLF) linebreaks in {0}".format(x)
		exit()
	
# parses command line arguments
def get_commandline_arguments ():
	
	import os
	import argparse
	
	parser = argparse.ArgumentParser()
	
	parser.add_argument("--z_option", required=True, dest="z_option", type=extant_file, help="name and path of the z_option file: 1st column barcode/sample name separated by tab from second column indicating the z_option: usually z0 or z1 or *blank*; if a sample missing fom this file it will be absent from newhybrids dataset", metavar="FILE")
	parser.add_argument("--vcf", required=True, dest="vcf", type=extant_file, help="vcf genotypes", metavar="FILE")

	args = parser.parse_args()
	print args
	
#	linebreak_check(args.popmapfile)
		
	return args

########


def check_congruence (z_optionfile, vcf):
	
	popmapsamples = []
	with open(z_optionfile, "r") as INFILE:
		for line in INFILE:
			if len(line) > 3:
				fields = line.strip("\n").split("\t")
				popmapsamples.append(fields[0])
	popmapsamples = set(popmapsamples)
	
	with open(vcf, "r") as INFILE:
		for line in INFILE:
			if line.startswith("##"):
				continue
			if line.startswith("#"):
				header_line = line.lstrip("#").strip("\n").split("\t")	
				break
	
	vcf_samples = header_line[9:]
	
	for samp in popmapsamples:
		if samp not in vcf_samples:
			print "vcf does not contain all z_optionfile samples"
			exit()
		
	print "all samples in z_optionfile are in .vcf, good to go!"


#######

def read_z_optionfile (z_optionfile):
	
	z_options = {}
	with open(z_optionfile,"r") as INFILE:
		for line in INFILE:
			if len(line) > 3:
				fields = line.strip("\n").split("\t")
				try:
					z_options[fields[0]] = fields[1]
				except IndexError:
					z_options[fields[0]] = ""
	
	return z_options



def vcf_to_newhybrids (vcf):
		
	with open(vcf, "r") as INFILE:
		for line in INFILE:
			if line.startswith("##"):
				continue
			if line.startswith("#"):
				header_line = line.lstrip("#").strip("\n").split("\t")
				break
	
	samples = header_line[9:]
	
	#	E3379_L96	43	.	ATCG	ATTG,ATCA,GTCA,GTCG	4179.3	.	AB=BLABLA	0/0:1:
	
	print samples
	genotype_dict = {x:[] for x in samples}
	print genotype_dict
	
	with open(vcf, "r") as INFILE:
		
		cnt = 0
		for line in INFILE:
			if line.startswith("#"):
				continue

			cnt += 1
			print str(cnt)

			fields = line.strip("\n").split("\t")

			for idx, sample in enumerate(samples):
				gt = fields[9+idx].split(":")[0].split("/")
				if not gt[0] == ".":
					gt_out = "".join([str(int(x) + 1) for x in gt])
				else:
					gt_out = "0"
				genotype_dict[sample].append(gt_out)	
				
			## first check that site occurs in the outgroup and fixed in the outgroup:
#			a = [ fields[x].split(":")[0].split("/") for x in popdict_vcf_idxes[outgr] ]
#			b = [x for genotype in a for x in genotype if x != "."]
#			pop_n = float( len( popdict_vcf_idxes[outgr] ) )*2
#			if len(b) >= presence_treshold_outgr * pop_n:
#				if len(set(b)) == 1:
					# we have a site that survives max_missing AND is fixed in the outgroup
					# now check the other populations:
#			print fields
	return genotype_dict

def filter_genotype_dict (genotype_dict, z_options):
	
	for k in genotype_dict.keys():
		if not k in z_options.keys():
			del genotype_dict[k]
	
	return genotype_dict
	
	

def write_newhybrids_format ( genotype_dict, output_file, z_options ):

 	line1 = "NumIndivs" + "\t" + str(len(genotype_dict.keys()))
	line2 = "NumLoci" + "\t" + str(len(genotype_dict.values()[0]))
	line3 = "Digits	1"
	line4 = "Format	Lumped"
 	line5 = "LocusNames" + "\t" + "\t".join( ["locus_"+str(x) for x in range(0,len(genotype_dict.values()[0])) ] )
	
	output_lines = [ line1, line2, line3, line4, line5 ]
	samples_map = []
	for idx, sample in enumerate(genotype_dict.keys()):
		outl = str(idx+1) + "\t" + z_options[sample] + "\t" + "\t".join( genotype_dict[sample] )
		samples_map.append( str(idx+1) + "\t" + sample  )
		output_lines.append(outl)
		
	with open(output_file, "w") as OUTF:
		OUTF.write( "\n".join(output_lines) + "\n"  )
	
	with open(output_file + ".samplemap.txt", "w") as OUTF:
		OUTF.write( "\n".join(samples_map) + "\n"  )

	with open(output_file + ".individuals.txt", "w") as OUTF:
		OUTF.write( "\n".join([x.split("\t")[1] for x in samples_map]) + "\n"  )

		
######################## MAIN


args = get_commandline_arguments ()
	
check_congruence (args.z_option, args.vcf)
	
z_options = read_z_optionfile (args.z_option) 
outfilename = args.vcf + ".newhybrids.txt"
datadict = vcf_to_newhybrids (args.vcf)
datadict = filter_genotype_dict (datadict, z_options)
write_newhybrids_format ( datadict, outfilename, z_options)

print "Done!"


