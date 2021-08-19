# Mathias Scharmann
# 
# 
# python 2.7




import argparse
import re
import subprocess
import os
import multiprocessing as mp

# diese Funktion gibt command line arguments mit flag an das script weiter
def get_commandline_arguments ():
	
	parser = argparse.ArgumentParser()
	
	parser.add_argument("-out", required=True,
	    dest="out", help="name of output file", metavar="FILE")
	
	args = parser.parse_args()
	
	return args


#####

def split_fasta_from_concatmatrix (taxon, genemodelfile, concat_seq):

	outlines = []
	with open(genemodelfile, "r") as INFILE:
		for line in INFILE:
			gene = line.split("=")[0]
			gene_start = int(line.strip("\n").split("=")[1].split("-")[0])-1 # for 0-based offset in python
			gene_end = int(line.strip("\n").split("=")[1].split("-")[1]) # if slicing in python [0:x], the value of index x will not be retrieved, but only everything before x
			
			#if [start:end], will retrieve items inclduding start until end-1
			
			outlines.append(">" + gene)
			outlines.append( concat_seq[gene_start:gene_end].replace("-", "") )	
			
			
			model_length = gene_end - (int(line.strip("\n").split("=")[1].split("-")[0]) -1)
			
		#	print "below_new"
		#	print model_length
		#	print len(concat_seq[gene_start:gene_end])
			
			#print ">" + gene 
			#print concat_seq[gene_start:gene_end].replace("-", "")
			
	with open("./tmp_working_dir/{0}.per_gene.cds.fasta".format(taxon), "w") as OUTFILE:
		OUTFILE.write("\n".join(outlines))

def read_phylip(INFILE):
	
	
	
	indict = {}
	with open(INFILE, "r") as infile:
		infile.readline() # remove header
		for line in infile:
			f = line.strip("\n")
			delim = re.compile('\s*')
			fields = delim.split(f)

			indict[fields[0]] = fields[1]
	
	print "read phylip input matrix"
	
	return indict

		
def read_all_genewise_fastas (all_genes, all_taxa):

	taxdict = {}
	for taxon in all_taxa:
		taxdict[taxon] = {}
		with open("tmp_working_dir/{0}.per_gene.cds.fasta".format(taxon), "r") as INFILE:
			for line in INFILE:
				if line.startswith(">"):
					ID = line.strip(">").strip("\n")
				else:
					seq = line.strip("\n")
					taxdict[taxon][ID] = seq

	# now, turn dict around to get per-locus data:
	genes_taxa_dict = {}
	for gene in all_genes:
		
		print gene
		
		genes_taxa_dict[gene] = {} # initialise second level
		
	
		for taxon in all_taxa:
			try:
				genes_taxa_dict[gene][taxon] = taxdict[taxon][gene]
			except KeyError:
				genes_taxa_dict[gene][taxon] = ""
		
	return genes_taxa_dict

##########

def MT_wrapper_muscle(cluster_dict, nthreads):

	print "there are {0} clusters".format(len(cluster_dict.keys()))
	
	results = {}
	
	
	pool = mp.Pool(nthreads) #use all available cores, otherwise specify the number you want as an argument
	for clust in sorted(cluster_dict.keys()):
		results[clust] = pool.apply_async(call_muscle_MT, args=(cluster_dict[clust], clust))
	pool.close()
	pool.join()
	
	# Get process results from the output queue
	#print output
	aligned_clusters = {}
	for idx, result in results.items():
#		print idx, result
		aligned_clusters[idx] = result.get()
#		print result.get()
	
	return aligned_clusters

def call_muscle_MT (seqs_dict, clust):
	
	print "\rrunning muscle on cluster \t", clust,
	
	alignment = {}
		
	with open("./tmp_working_dir/"+ "{0}.for_muscle.fasta".format(clust), "w") as OUTFILE:
		for k, v in seqs_dict.items():
			OUTFILE.write(">"+str(k)+"\n"+v+"\n")
	
	bash_command = "muscle -quiet -in ./tmp_working_dir/" + "{0}.for_muscle.fasta".format(clust) # -quiet
#	bash_command = "mafft-linsi --quiet {0}.for_muscle.fasta".format(clust)

	p = subprocess.Popen(bash_command, shell=True, stdout=subprocess.PIPE) #, stderr=subprocess.STDOUT)
	
	rec = "_"
	
	while True:
		line = p.stdout.readline()
#		print line,
		if line == '' and p.poll() != None:
			# deposit also last line of file:
			alignment[taxon] = [rec]
			# after thatt break
			break
		if line.startswith(">"): # next taxon
			
			# first round:
			if rec == "_":
				taxon = line.strip(">").strip("\n")
				rec = ""
			
			if rec != "":
#				print "deposting", taxon
#				print rec
				# make deposit of previous rec, reset rec, adopt new taxon
				alignment[taxon] = [rec]
				# after deposit, reset rec to empty		
				rec = ""
				taxon = line.strip(">").strip("\n")
				
			# reset taxon in case record is empty
			else:
				taxon = line.strip(">").strip("\n")
				rec = ""
		else: # same taxon, record sequence
			rec += line.strip("\n")	
		
#	print alignment[cluster]
	
	os.remove("./tmp_working_dir/{0}.for_muscle.fasta".format(clust))
						
	out_dict = {}
	for taxon, seqlist in alignment.items():
		out_dict[taxon] = "".join(seqlist)				
		
	return out_dict
	
############

try:
	os.mkdir("./tmp_working_dir")
except OSError:
	None
	
print "working in tmp_working_dir"

phymatrix = read_phylip("ortho_MI.supermatrix.4.phy")
genemodelfile = "ortho_MI.supermatrix.4.model"

#phymatrix = read_phylip("fake.phy")
#genemodelfile = "fake.model.txt"

for taxon in phymatrix.keys():
	split_fasta_from_concatmatrix (taxon, genemodelfile, phymatrix[taxon])

print "split concatenated matrix per gene, written to file"

# now we have for each taxon a fasta file. -- the one to be added comes in now:
all_genes = set()
with open(genemodelfile, "r") as INFILE:
	for line in INFILE:
		all_genes.add( line.split("=")[0] )

#all_genes = set(["WAG,cc4779-1.inclade1.ortho1.aln-cln_gene1"])#, 
#all_genes = list(all_genes)[:100]

all_taxa = phymatrix.keys()

all_taxa = ["spis",
"lini",
"muir",
"xant",
"brun",
"dubi",
"rubr",
"olen",
"plat",
"eric",
"Leuc",
"Leucadendron_argenteum",
"Leucospermum_erubescens",
"Leucadendron_ericifolium",
"Mimetes_cucullatus",
"Leucadendron_galpinii"
]

#all_taxa = ["Dionaea", "khasiana","WQUF", "Brap","CUTE","FUENF"]

# a dictionary with 2 levels {gene:{taxon:seq}}
genes_taxa_dict = read_all_genewise_fastas (all_genes, all_taxa)

# until here data are correct 

nthreads = 12
aligned_genes = MT_wrapper_muscle(genes_taxa_dict, nthreads)



print "\nmuscle done"

# wrap up

actually_used_genes = []
for gene in sorted(all_genes):
	if len(aligned_genes[gene].keys()) > 0:
		actually_used_genes.append(gene)

new_supermatrix = {}
for taxon in all_taxa:
	new_supermatrix[taxon] = ""
	for gene in actually_used_genes:
		try:
			seq = aligned_genes[gene][taxon]
		except KeyError:
#			print aligned_genes[gene].keys()
			seq = "-"*len( aligned_genes[gene][ aligned_genes[gene].keys()[0] ] )
		new_supermatrix[taxon] += seq
							

# columns to rows/lines:
outlines = []
	
ntaxa = len(all_taxa)
len_align = len( new_supermatrix[new_supermatrix.keys()[0]] )
header = str(ntaxa) + " " + str(len_align)
outlines.append(header)

for taxon in sorted(all_taxa):
	outlines.append( taxon + "    " + new_supermatrix[taxon] )
	
with open("new_supermatrix.phy", "w") as OUTFILE:
	OUTFILE.write("\n".join(outlines))


# build the RAxML gene model file (-q option):
total_length = 0
outlines = []
for gene in actually_used_genes:
	genelength = len( aligned_genes[gene][ aligned_genes[gene].keys()[0] ] )
	firstpos = total_length + 1
	lastpos = total_length + genelength
	outlines.append( gene + "=" + str(firstpos) + "-" + str(lastpos) )
	total_length += genelength
with open("new_supermatix.model", "w") as OUTFILE:
	OUTFILE.write("\n".join(outlines))




print "Done!"

	
