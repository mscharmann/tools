
import sys

# DNA codon table from BioPython 
# https://github.com/biopython/biopython/blob/master/Bio/Data/CodonTable.py
# Data from NCBI genetic code table version 4.2

codon_table={
'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
'TAT': 'Y', 'TAC': 'Y', 'TGT': 'C', 'TGC': 'C',
'TGG': 'W', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L',
'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q',
'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R',
'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T',
'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K',
'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R',
'AGG': 'R', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V',
'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A',
'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E',
'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
'GGG': 'G', 'TAA': '*', 'TAG': '*', 'TGA': '*'}



def read_fasta (infile):

	fasta_dict = {}
	seq = ""
	name = "dummy"
	with open(infile, "r") as F:
		for line in F:
			if line.startswith(">"):
				fasta_dict[name] = seq
				name = line.lstrip(">").rstrip("\n")
				seq = ""
			else:
				seq += line.strip("\n")

	# last record:
	fasta_dict[name] = seq
	del fasta_dict["dummy"]
	return fasta_dict


def write_fasta (seqdict, outfile):

	outlines = []
	the_order = sorted( seqdict.keys() )
	for id in the_order:
		outlines.append(">" + id)
		outlines.append(seqdict[id])
	
	with open(outfile, "w") as F:
		F.write("\n".join(outlines) + "\n")



def make_unique_pairs (inlist):
	
	pairs = set()
	for a in inlist:
		for b in inlist:
			if a!= b:
				pair = "XjXjX".join( sorted([a,b]) )
				pairs.add( pair )
	return pairs



def simpleclean(infile, window_size, window_threshold):
	
	seq_dict = read_fasta (infile)
	
	all_pairs = make_unique_pairs (seq_dict.keys())

	pairwise_dists_per_column = {} # counting non-synonymous changes per codon 
	for pair in all_pairs:
		pairwise_dists_per_column[pair] = []
		[a,b] = pair.split("XjXjX")
		seq_a = seq_dict[a].upper()
		seq_b = seq_dict[b].upper()
	#	print seq_a, seq_b
		idx = 0
		while idx < len( seq_a ):
			ca = "".join([ x for x in seq_a[idx:(idx+3)] if x not in ["-","N","n","?"]])
			cb = "".join([ x for x in seq_b[idx:(idx+3)] if x not in ["-","N","n","?"]])
			if len(ca) == 3 and len(cb) == 3:
				if codon_table[ca] == codon_table[cb]:
					pairwise_dists_per_column[pair].append(0)
				else:
					pairwise_dists_per_column[pair].append(1)
			else:
				pairwise_dists_per_column[pair].append(0)
			idx += 3
		

	id_sums_per_column = {} # average pairwise distances (all equal / Jukes-Cantor) ; missing data / gap counted as identical
	for id in seq_dict.keys():
		id_sums_per_column[id] = [0]*len( pairwise_dists_per_column[pairwise_dists_per_column.keys()[0]] ) 
		ids_pairs = [ x for x in all_pairs if id in x ]
		for i in range(len( pairwise_dists_per_column[pairwise_dists_per_column.keys()[0]] )): # loop over columns
			for p in ids_pairs:
				id_sums_per_column[id][i] += pairwise_dists_per_column[p][i]
			id_sums_per_column[id][i] = float( id_sums_per_column[id][i] ) / float(len( ids_pairs ))	# divide by number of pairs	
	

	# now window-slide & mask:
	masked_seqs_dict = {}
	for id in id_sums_per_column.keys():
		to_be_masked = []
		idx = -1
		mismatches = id_sums_per_column[id]
		seq = seq_dict[id]
		for i in range(len( pairwise_dists_per_column[pairwise_dists_per_column.keys()[0]] )):
			if i < len( pairwise_dists_per_column[pairwise_dists_per_column.keys()[0]] ) - window_size:
				idx += 1
				start_idx = idx 
				end_idx = start_idx + window_size
				mismatch_window = mismatches[start_idx:end_idx]
				if sum(mismatch_window) >= window_threshold:
					to_be_masked.append( [start_idx*3, end_idx*3] )
		outseq = list(seq)
		for bad_window in to_be_masked:
			for i in range(bad_window[0],bad_window[1]):
				outseq[i] = "N"
		masked_seqs_dict[id] = "".join(outseq)
	
	return masked_seqs_dict



####

if len(sys.argv) != 4:
	print "usage: python simpleclean.py cds_aln_file window_size window_threshold"
	exit()


print "masking residues by average pairwise distance in sliding windows"	
result = simpleclean(sys.argv[1],int(sys.argv[2]),int(sys.argv[3]))
write_fasta (result, sys.argv[1] + ".masked.aln")
print "Done!"


