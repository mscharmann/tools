# Python 2.7.6
# fasta_to_majority_consensus.py
# Mathias Scharmann

# 3 July 2015

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
	
	
	outseq = []
	if len(seqs) > 0:
		for column_index in range(0, len(seqs[0])): # to loop over all columns
			aligned_column = "".join([ seq[column_index] for seq in seqs ])
#			print aligned_column
			maxcount = 0
			for nuc in ["A","T","G","C","-"]:
				if aligned_column.count(nuc) > maxcount:
					maxcount = aligned_column.count(nuc)
					consensus = nuc
			
			outseq.append(consensus)
		outseqstring = "".join(outseq)

		header_line = ">{0}\n".format(name)
		seq_line = outseqstring + "\n"
		FASTA.write(header_line)
		FASTA.write(seq_line)	
       	

print "Done!"
	
	