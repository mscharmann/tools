#######################################
"""
aim: for design of exon-exon junction primers, spanning the junction, this website offers a nice solution:


1) get spliced version of gene with exon coordinates projected onto it with gffread -x -W
			
			
2) open the resulting fasta, and mark the exon-exon junctions IN THE CODING SEQUENCE with a "-" character (keep track of the offset!)
3) got to primer 3 website which has option "SEQUENCE_OVERLAP_JUNCTION_LIST"
         https://primer3.ut.ee/

example:

python mark_exon-exon_junctions_in_CDS.py spliced_with_coords.fa

"""




import sys



def read_wrapped_or_unwrapped_fasta (infile):
	outlines = []
	with open(infile, "r") as INFILE:
		first_id = INFILE.readline().strip("\n").strip(">")
		outlines.append(first_id)
		seq = ""
		for line in INFILE:
			line_clean = line.strip("\n")
			if line_clean.startswith(">"):
				outlines.append(seq)
				outlines.append(line_clean.strip(">"))
				seq = ""
			else:
				if len(line_clean) > 0:
					seq += line_clean
		# append last seq
		outlines.append(seq)
	i=0
	j=1
	out_dict = {}
	for x in range(len(outlines)/2):
		out_dict[outlines[i]] = outlines[j]
		i += 2
		j += 2
	return out_dict


	


######


infile = sys.argv[1]

inseqs = read_wrapped_or_unwrapped_fasta(infile)

###

outlines = []
for name, seq in inseqs.items():
	seq_with_marked_exons = ""
	print(name, len(seq))
	exons1 = name.split("segs:")[1].split(",")
	exons2 = [x.split("-") for x in exons1]
	print (exons2)  ## 1-based and inclusive
	for x in exons2[:-1]:
		print(x)
		exons_seq = seq[int(x[0])-1:int(x[1])]
		print( len(exons_seq), int(x[1])-int(x[0]) )
		seq_with_marked_exons += exons_seq
		seq_with_marked_exons += "-"
	# add the last one
	exons_seq = seq[int(exons2[-1][0])-1:int(exons2[-1][1])]
	seq_with_marked_exons += exons_seq
	print(seq_with_marked_exons)
	outlines.append(">" + name)
	outlines.append(seq_with_marked_exons)
	

with open(infile + ".exon_junctions_marked_in_sequence.fasta", "w") as O:
	O.write("\n".join(outlines) + "\n")
	
	
