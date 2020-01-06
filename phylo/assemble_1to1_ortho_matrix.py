## python 3

## python assemble_1to1_ortho_matrix.py PATH FILE_SUFFIX
# reads fasta alignments per-gene from a directory and assemble a supermatrix
# please check that fasta headers conform to >taxon_blabla
# name of the fasta file = name of gene in supermatrix model file

import sys, os

allfiles = os.listdir(sys.argv[1])
infiles = sorted([x for x in allfiles if x.endswith(sys.argv[2])])

print("number of input files:")
print(len(infiles))
print("parsing files to get all taxa")

all_taxa = set()
for f in infiles[:1000]:
	with open(sys.argv[1] + f, "r") as I:
		for line in I:
			if line.startswith(">"):
				t = line.strip(">").split("_")[0]
				all_taxa.add(t)

all_taxa = list(sorted(all_taxa))				
print("taxa:")
print (all_taxa)
print("second pass, now to get sequences and assembling outputs")


one_to_one_count = 0
out_seq_dict = {t:[] for t in all_taxa}
modelfile_lines = []
left_idx = 1
right_idx = 0

for f in infiles:
	nseq = 0
	taxa = []
	with open(sys.argv[1] + f, "r") as I:
		for line in I:
			if line.startswith(">"):
				nseq += 1
				t = line.strip(">").split("_")[0]
				taxa.append(t)
#	print(f, nseq, taxa)
	if list(sorted(taxa)) == all_taxa: # taxa are complete
		if nseq == len(taxa): # unique sequence per tax
			print(f)
			one_to_one_count += 1 
			# get the seqs
			with open(sys.argv[1] + f, "r") as I:
				for line in I:
					if line.startswith(">"):
						t = line.strip(">").split("_")[0]
					else:
						seq = line.strip("\n")
						out_seq_dict[t].append( seq )
				
			# prep modelfile
			aligned_len = len(seq)
			right_idx += aligned_len
			mline = "DNA" + "," + f.strip(sys.argv[2]) + "=" + str( left_idx ) + "-" + str( right_idx )
			left_idx += aligned_len
			modelfile_lines.append( mline )
				

with open("supermatrix.model","w") as outfile2:
		outfile2.write( "\n".join( modelfile_lines ) + "\n")


for k,v in out_seq_dict.items():
	print(k,len(v),len("".join(v))) 

outlines_fasta = [">" + k + "\n" + "".join(v) for k,v in out_seq_dict.items()]
with open("supermatrix.fasta","w") as outfile1:
		outfile1.write( "\n".join( outlines_fasta ) + "\n")


print("exported " + str(one_to_one_count) + " 1-to-1 orthologous genes, DONE!")
