#!/usr/local/bin/python
# Python 3

# Mathias Scharmann




#import scipy.special
import sys

lower_threshold = float( sys.argv[1] )
upper_threshold = float( sys.argv[2] )


################################## MAIN

for line in sys.stdin:
	if line.startswith("#"):
		sys.stdout.write(line)
		continue
	if len(line) < 2: # empty lines or so
		sys.stdout.write(line)
		continue
	fields = line.strip("\n").split("\t")
	gt_fields = fields[9:]
	gts = "".join( [x.split(":")[0] for x in gt_fields] ).replace("/","")		
	try:
		alleles = set(gts)
		alleles.discard(".")
		n_alleles = len(alleles)
	except TypeError:
		n_alleles = 3 # dummy

	if n_alleles == 1: # fixed sites pass
		sys.stdout.write(line)
		continue
				
	elif n_alleles == 2:	# only bi-allelic sites, NOT MISSING SITES
		# get heterozygote genotypes
		gts_all = [x.split(":")[0].split("/") for x in gt_fields]
		gts_het_idxes = [idx for idx,x in enumerate(gts_all) if x[0] != x[1]]
		
		#gts_het_idxes = [gts_all.index(x) for x in gts_all if x[0] != x[1]]
		if len(gts_het_idxes) < 1: # sites whithout any heterozygotes pass
			sys.stdout.write(line)
			continue
		else: # get allelic depths in heterozygotes
			#print("HERE1", line)
			#print(gts_het_idxes)
			allelic_depths = [gt_fields[i].split(":")[3].split(",") for i in gts_het_idxes]
			#print(line, gts_all, allelic_depths)
			count_1 = 0
			count_2 = 0
			for d in allelic_depths:
				count_1 += int(d[0])
				count_2 += int(d[1])
			AB = float(count_1)/(count_1+count_2)
			if lower_threshold < AB < upper_threshold:
				sys.stdout.write(line) # sites with heterozygotes only pass if AB within these bounds



