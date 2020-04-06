# common_loci_finder.py list_of_infiles common_to_howmany
import itertools
import sys


file_list_of_infiles = sys.argv[1]
common_to_howmany = int(sys.argv[2])

with open( file_list_of_infiles, "r" ) as INFILE:
	infile_list = [ l.strip("\n") for l in INFILE ]

	
data_dict = {}
for infile in infile_list:
	with open( infile, "r" ) as INFILE:
		data_dict[infile] = set( [ l.strip("\n") for l in INFILE ] )


import itertools
all_combs = list( itertools.combinations(data_dict.keys(), common_to_howmany)	)	


# now find the set of loci that are shared in each combination:
master_set = set()
for comb in all_combs:
	
	combs_set = data_dict[comb[0]]
	for i in comb[1:]:
		combs_set = combs_set.intersection( data_dict[i] )
	master_set = master_set.union( combs_set )


with open("common_loci", "w") as OUTFILE:
	OUTFILE.write("\n".join( sorted(list(master_set)) ))