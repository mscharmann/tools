# readfq_filter_phred20_in_barcode.PE_reads.py

# Mathias Scharmann
# fq parser downloaded
# 08 June 2018
# python 2.7

# 
# python ../tools/readfq_filter_phred20_in_barcode.PE_reads.py try_R1.fq.gz try_R2.fq.gz

n_leading_bases_to_check = 8

import sys, gzip


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

	first_reads = sys.argv[1]
	second_reads = sys.argv[2]
	
	outfile_firstreads = first_reads + ".HQbarcode.fq"
	outfile_secondreads = second_reads + ".HQbarcode.fq"
	
	n_survived, n_dropped = 0, 0
	bad_entries = []
	total_entry_count = 0
	output_buffer = []
#	with gzip.open(outfile_firstreads, "w") as OUTF:
#		with gzip.open(first_reads, "r") as INF:
	with open(outfile_firstreads, "w") as OUTF:
		with open(first_reads, "r") as INF:

			for name, seq, qual in readfq(INF):
				total_entry_count += 1
				qual_numeric = [ ord(x) - 33 for x in qual.rstrip('\n') ] ## from http://pythonforbiologists.com/index.php/business-card-fastq-parser-exercise/
				if len( [ x for x in qual_numeric[:n_leading_bases_to_check] if x >= 20] ) == n_leading_bases_to_check:
					n_survived += 1
					outentry = "@" + name + "\n" + seq + "\n+\n" + qual + "\n"
					output_buffer.append( outentry )
#					print total_entry_count, n_survived
				else:
					n_dropped += 1
					bad_entries.append( total_entry_count )
				if (len(output_buffer)/100000.0).is_integer():
					OUTF.write( "".join(output_buffer) )
					print "output written"
					output_buffer = []
		OUTF.write( "".join(output_buffer) )
		print "last written R1, now filtering R2"
		output_buffer = []
	
	
	bad_entries_set = set( bad_entries )
	total_entry_count = 0
	output_buffer = []
#	with gzip.open(outfile_secondreads, "w") as OUTF:
#		with gzip.open(second_reads, "r") as INF:
	with open(outfile_secondreads, "w") as OUTF:
		with open(second_reads, "r") as INF:
			for name, seq, qual in readfq(INF):
				total_entry_count += 1
				if total_entry_count not in bad_entries_set:
					outentry = "@" + name + "\n" + seq + "\n+\n" + qual + "\n"
					output_buffer.append( outentry )
				if (len(output_buffer)/100000.0).is_integer():
					OUTF.write( "".join(output_buffer) )
					print "output written"
					output_buffer = []
		OUTF.write( "".join(output_buffer) )
		print "last written R2, now done"
		output_buffer = []
					
		
	sys.stderr.write( "total:		" + str(n_dropped + n_survived) + "\n" + "passed:		" + str(n_survived) + "\n" + "dropped:	" + str(n_dropped) + "\n" + "prop dropped:	" + str( float(n_dropped)/float( n_dropped + n_survived ) ) + "\n" )
				

	

