# Python 2.7.6
# fasta_to_majority_consensus.py
# Mathias Scharmann

# 14 June 2016

"""
example usage:
cat blasthits_557475_L96.fasta | python $tb/fasta_to_majority_consensus.py -out fufufu.txt

"""


import argparse
import sys

#######


def get_commandline_arguments ():
	
	parser = argparse.ArgumentParser()
		
	parser.add_argument("-out", required=True,
	    dest="out", help="name of output file", metavar="FILE")
		
	args = parser.parse_args()
	
	return args

revcombs = {'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'K': 'GT', 'M': 'AC', 'N': 'ACGT', 'S': 'CG', 'R': 'AG', 'W': 'AT', 'V': 'ACG', 'Y': 'CT', 'A' : 'A', 'C':'C', 'T':'T', 'G':'G', 'n':'-', '-':'-'}


def get_IUPAC_amb_code (nucs):
	
	# nucs must be a "".joined string of the sorted list of nucleotides
	
	# the IUPAC ambiguity codes:
	# the keys in this IUPAC ambiguity code dictionary are produced by:
	# "".join(sorted(list("NUCLEOTIDES"))) -> consistent lookup possible without additional loops!
	combs = {'AC':'M', 'GT':'K', 'CG':'S', 'AT':'W', 'AG':'R', 'CT':'Y',
	'ACG':'V', 'CGT':'B', 'AGT':'D', 'ACT':'H', 'ACGT':'N',
	'A' : 'A', 'C':'C', 'T':'T', 'G':'G', 'n':'-', '-':'-', 'N':'N'}

	
	out_code = combs[nucs]
	
	return out_code
	
def make_cluster_consensus(cluster_dict):
	
	# make consensus sequence out of individual sequences in the binned_dict:
	consensus_dict = {}
	for cluster, seqs in cluster_dict.items():
		outseq = []
		print "working on cluster {0}".format(cluster)
#		print seqs[0]
		if len(seqs) > 0:
			for column_index in range(0, len(seqs[0])): # to loop over all columns
				aligned_column = "".join([ seq[column_index] for seq in seqs ])
	#			print aligned_column
				maxcount = 0
				for nuc in ["A","T","G","C"]:
					if aligned_column.count(nuc) > maxcount:
						maxcount = aligned_column.count(nuc)
						consensus = nuc
				
			
				outseq.append(consensus)
		outseqstring = "".join(outseq)
#		print outseqstring
		consensus_dict[cluster] = outseqstring

	return consensus_dict


#####

def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break

def make_amb_seqs (indict):
		
	# make consensus sequence out of individual sequences in the binned_dict:
	seqs = indict.values()
		
	outseq = []
	for column_index in range(0, len(seqs[0])): # to loop over all columns
#		print column_index
		aligned_column = [ seq[column_index] for seq in seqs ]
		resolved = [ revcombs[x] for x in aligned_column ] 
		aligned_column = list("".join(resolved))
		nucs = [x for x in aligned_column] # remove all the n -
		if "-" not in nucs: # retain gaps as they are! Do not fill with seq from other fasta entries. consensus of a gap is always a gap.
			if len(nucs) == 0:
				outnuc = "-"
			else:
				outnuc = get_IUPAC_amb_code( "".join(sorted(list(set(nucs)))) )
		else:
			outnuc = "-"
		print column_index, "\r",	
		outseq.append(outnuc)
	outseqstring = "".join(outseq)

	return outseqstring


	
######### MAIN ########


	
if sys.stdin.isatty(): # returns False if something is in STDIN
	print "no input data piped, dying"
	exit()

args = get_commandline_arguments ()

with open(args.out, "w") as FASTA:
	n, slen, qlen = 0, 0, 0
	
	fastadict = {}
	seqs = []
	for name, seq, qual in readfq(sys.stdin):
		fastadict[name] = seq
		seqs.append(seq)
	
	outseqstring = make_amb_seqs (fastadict)

	header_line = ">ambiguity_consensus\n"
	seq_line = outseqstring + "\n"
	FASTA.write(header_line)
	FASTA.write(seq_line)	
       	

print "Done!"
	
	