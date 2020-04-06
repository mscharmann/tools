"""
Input:
- aligment in phylip format
- raxml partition format file

Output: 
- a directory with trees estimated separately for each partition!


module load java
module load gdc
module load mafft python/2.7
PATH=$PATH:/cluster/project/gdc/people/schamath/tools

"""

import sys, os
import subprocess

###

def read_phylip(INFILE):
	
	indict = {}
	with open(INFILE, "r") as infile:
		infile.readline() # remove header
		for line in infile:
			fields = line.strip("\n").split( )
			indict[fields[0]] = fields[1]

	return indict


def write_phylip (concat_dict, outfile):
	
	with open(outfile, "w") as OUTFILE:
		ntaxa = len(concat_dict.keys())
		len_align = len(concat_dict[concat_dict.keys()[0]]) # gets value of first element in dictionary -> this OK since all seqs have same length
		header = str(ntaxa) + " " + str(len_align) + "\n"
		OUTFILE.write(header)
		for sample, seq in concat_dict.items():
			out_line = sample + "    " + seq + "\n"
			OUTFILE.write(out_line)
	OUTFILE.close()


def split_alignment ( alignment_file, partitions_file, outdir ):
	
	master_aln = read_phylip( alignment_file )
	
	partitions_dict = {}
	with open(partitions_file, "r") as INFILE:
		for line in INFILE:
			name = line.split(",")[1].split("=")[0]	
			start_idx = int(line.split("=")[1].split("-")[0]) -1 # zero-based offset for python idxes
			end_idx = int(line.split("=")[1].split("-")[1] ) + 1 # string/list slicing in python: right end idx NOT included in slice
			partitions_dict[name] = [start_idx, end_idx]								
	
	cnt = 0
	for name, idxes in partitions_dict.items():
		cnt += 1		
		outdict = {tax: master_aln[tax][partitions_dict[name][0] : partitions_dict[name][1]] for tax in master_aln.keys() }		
		write_phylip (outdict , outdir + "/" +  name + ".phy")
		
		raxml_run_name = "run_" + str(cnt) + "_" 
		call_raxml (4, raxml_run_name, name + ".phy", outdir )
		
		os.remove( outdir + "/" +  name + ".phy" )
		os.remove( outdir + "/" + "RAxML_parsimonyTree." + raxml_run_name )
		os.remove( outdir + "/" + "RAxML_result." + raxml_run_name )
		os.remove( outdir + "/" + "RAxML_log." + raxml_run_name )
		os.remove( outdir + "/" + "RAxML_info." + raxml_run_name )
		try:
			os.remove( outdir + "/" +  name + ".phy.reduced"  )
		except OSError:
			None
		os.rename( outdir + "/" + "RAxML_bestTree." + raxml_run_name, outdir + "/" + name + ".raxml.tre" )

def call_raxml (num_cores, raxml_run_name, aln_file, outdir):
	
	cmd = " ".join ( ["cd", outdir,";", "raxml","-T",str(num_cores),"-p","12345","-s",\
			   aln_file,"-n", raxml_run_name ,"-m GTRCAT"] )
	print cmd
	return_code = subprocess.call(cmd, shell=True)  
#	p = subprocess.Popen(cmd,stdout=subprocess.PIPE)
#	out = p.communicate()



########## main

if len(sys.argv) != 4:
	print "usage: python raxml_for_each_partition.py alignment_file partitions_file outdir"
	exit()

alignment_file = sys.argv[1]
partitions_file = sys.argv[2]
outdir = sys.argv[3]
try:
	os.mkdir( outdir )
except OSError:
	None

split_alignment( alignment_file, partitions_file, outdir )

