# readfq_filter_phred20_in_barcode.py

# Mathias Scharmann
# fq parser downloaded
# 04 Oct 2016
# python 2.7

# pipe *.fq files into script, filters for phred >= 20 in each of the first 5 bases, writes surviving entries to STDOUT
# example usage :
#	zcat 20151211.A-ddRAD_Nph_1_R1.adaptertrim.fq.gz | head -1000000 | python readfq_filter_phred20_in_barcode.py | gzip > fufu.fq.gz



import sys


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


############

if __name__ == "__main__":
	
	if sys.stdin.isatty(): # returns False if something is in STDIN
		print "no input data piped, dying"
		exit()
	
	n_survived, n_dropped = 0, 0
	seen_seqs = set()
	seqs_counts = {}
	for name, seq, qual in readfq(sys.stdin):
		qual_numeric = [ ord(x) - 33 for x in qual.rstrip('\n') ] ## from http://pythonforbiologists.com/index.php/business-card-fastq-parser-exercise/
		if len( [ x for x in qual_numeric[:5] if x >= 20] ) == 5:
			n_survived += 1
			outentry = "@" + name + "\n" + seq + "\n+\n" + qual + "\n"
			sys.stdout.write( outentry )
			sys.stdout.flush()
		else:
			n_dropped += 1
			
	
	sys.stderr.write( "total:		" + str(n_dropped + n_survived) + "\n" + "passed:		" + str(n_survived) + "\n" + "dropped:	" + str(n_dropped) + "\n" + "prop dropped:	" + str( float(n_dropped)/float( n_dropped + n_survived ) ) + "\n" )
				

	

