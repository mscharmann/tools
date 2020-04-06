# raxml_for_each_partition_multithreading.py
# Mathias Scharmann
# 2018-02-21



"""
Input:
- aligment in phylip format
- raxml partition format file

Output: 
- a directory with trees estimated separately for each partition!


module load gdc
module load python/2.7
PATH=$PATH:/cluster/project/gdc/people/schamath/tools


"""

import sys, os
import subprocess

###

def is_non_zero_file(fpath):  
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0

def read_phylip(INFILE):
	
	indict = {}
	with open(INFILE, "r") as infile:
		infile.readline() # remove header
		for line in infile:
			fields = line.strip("\n").split( )
			indict[fields[0]] = fields[1]

	return indict


def write_phylip_clean (cl_dict, outfile):
		
	with open(outfile, "w") as OUTFILE:
		ntaxa = len(cl_dict.keys())
		len_align = len(cl_dict[cl_dict.keys()[0]]) # gets value of first element in dictionary -> this OK since all seqs have same length
		header = str(ntaxa) + " " + str(len_align) + "\n"
		OUTFILE.write(header)
		for sample, seq in cl_dict.items():
			out_line = sample + "    " + seq + "\n"
			OUTFILE.write(out_line)
	OUTFILE.close()


def split_alignment_and_do_raxml ( alignment_file, partitions_file, outdir, parallel_processes ):
	
	master_aln = read_phylip( alignment_file )
	
	partitions_dict = {}
	with open(partitions_file, "r") as INFILE:
		for line in INFILE:
			name = line.split(",")[1].split("=")[0]	
			start_idx = int(line.split("=")[1].split("-")[0]) -1 # zero-based offset for python idxes
			end_idx = int(line.split("=")[1].split("-")[1] ) #+ 1 # string/list slicing in python: right end idx NOT included in slice
			partitions_dict[name] = [start_idx, end_idx]								
	
	cnt = 0
	processes = set()
	max_processes = parallel_processes
	num_cores_per_process = num_cores = 2
	for name, idxes in partitions_dict.items():
		cnt += 1		
		outdict = {tax: master_aln[tax][partitions_dict[name][0] : partitions_dict[name][1]] for tax in master_aln.keys() }		
		# drop taxa that have only undetermined characters "-", otherwise RAxML will complain
		clean_dict = {k:v for k,v in outdict.items() if set(v) != set("-") }
		
		# sometimes, alignment is slightly incorrect and contains ONLY missing data for a locus; catch these errors and continue. Also, gene trees aren't informative if there is not at least 4 tips on them.
		if len(clean_dict.keys()) >= 4:
			
			# also, don't run again if there is already a treefile with content in it:
			if not is_non_zero_file("./" + outdir + "/" + name + ".raxml.tre"):

				write_phylip_clean (clean_dict , outdir + "/" +  name + ".phy")
		
				raxml_run_name = "run_" + str(cnt) + "_"
				cmd = " ".join ( ["cd", outdir,";", "raxml","-T",str(num_cores),"-p","12345","-s",\
								name + ".phy" ,"-n", raxml_run_name ,"-m GTRCAT", "--silent" , ";",\
								"rm", name + ".phy", "RAxML_parsimonyTree." + raxml_run_name,\
								"RAxML_result." + raxml_run_name,\
								"RAxML_log." + raxml_run_name,\
								"RAxML_info." + raxml_run_name,\
								name + ".phy.reduced",\
								";", "mv", "RAxML_bestTree." + raxml_run_name, name + ".raxml.tre"] )
	
	#			print cmd
				processes.add(subprocess.Popen([cmd, raxml_run_name], shell = True))
				if len(processes) >= max_processes:
					os.wait()
					processes.difference_update([ p for p in processes if p.poll() is not None])

########## main

if len(sys.argv) != 5:
	print """usage: python raxml_for_each_partition_multithreading.py alignment_file partitions_file outdir total_cores\n
RAxML will be called on 2 cores for each parallel job"""
	
	exit()

alignment_file = sys.argv[1]
partitions_file = sys.argv[2]
outdir = sys.argv[3]
try:
	os.mkdir( outdir )
except OSError:
	None

total_cores = int( sys.argv[4] )
num_cores_per_process = 2
parallel_processes = total_cores / num_cores_per_process

split_alignment_and_do_raxml ( alignment_file, partitions_file, outdir, parallel_processes )

