# readfq_filter_phred20_in_barcode.PE_reads.py

# Mathias Scharmann
# fq parser downloaded
# 11 October 2018
# python 2.7

# 
# python readfq_filter_phred20_in_barcode.demultiplex.qual_filter.PE_reads.py --r1 lib_play.R1.fq --r2 lib_play.R2.fq --rmot1 TGCAGG --rmot2 TAA --w 0.15 --s 0 --b barcodes.txt --o samples




import sys, gzip, os
import argparse

# this function checks if file exists and exits script if not
def extant_file(x):
	"""
	'Type' for argparse - checks that file exists but does not open.
	"""
	
	if not os.path.exists(x):
		print "Error: {0} does not exist".format(x)
		exit()
	x = str(x)
	return x

# parses arguments from command line
def get_commandline_arguments ():
	
	parser = argparse.ArgumentParser()

	parser.add_argument("--r1", required=True, help="first read file", type=extant_file, metavar="FILE")
	parser.add_argument("--r2", required=True, help="second read file", type=extant_file, metavar="FILE")	
	parser.add_argument("--rmot1", required=True, help="restriction motif on first read", metavar="STRING")
	parser.add_argument("--rmot2", required=True, help="restriction motif on second read", metavar="STRING")
	parser.add_argument("--w", required=True, help="sliding window length, in prop. of read length", metavar="FLOAT")
	parser.add_argument("--s", required=True, help="min avg quality in window", metavar="FLOAT")
	parser.add_argument("--b", required=True, help="barcodes file", type=extant_file, metavar="FILE")
	parser.add_argument("--o", required=True, help="output directory", metavar="PATH")
	args = parser.parse_args()
	
	# finish
	return args


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


def read_barcodefile (infile):
	
	barcodes = []
	with open(infile, "r") as F:
		for line in F:
			if len(line) > 0:
				barcodes.append(line.strip("\n"))
	
	bc_dict = {}
	for b in barcodes:
		try:
			bc_dict[len(b)].append(b)
		except KeyError:
			bc_dict[len(b)] = [b]
	return bc_dict


############

if __name__ == "__main__":
	
	
	args = get_commandline_arguments ()
	
	first_reads = args.r1
	second_reads = args.r2
	
	barcodes_dict = read_barcodefile(args.b)
	try:
		os.mkdir(args.o)
	except OSError:
		print "outfile directory ", args.o, "already exists"
		None

	try:	
		os.remove(first_reads + ".tmp")
	except OSError:
		None
	try:
		os.symlink(first_reads, first_reads + ".tmp")
	except OSError:
		None				
	try:
		os.remove(second_reads + ".tmp")
	except OSError:
		None
	try:	
		os.symlink(second_reads, second_reads + ".tmp")
	except OSError:
		None
	try:
		os.remove(first_reads + ".discard")
	except OSError:
		None
	try:
		os.remove(second_reads + ".discard")
	except OSError:
		None
	
	logfile = open("demultiplex.log.txt", "a")
	
	with open(first_reads + ".tmp", "r") as INF1:
		INF1.readline()
		readlen = len(INF1.readline().strip("\n"))
	window_size = int( readlen*float( args.w ) )
	min_qual_sum = window_size*int(args.s)
	
	output_files_dict = {}		
	for k,v in barcodes_dict.items():
		for b in v:
			output_files_dict[b+"1"] = open(args.o + "/" + "sample_" + b + ".1.fq" , "a")
			output_files_dict[b+"2"] = open(args.o + "/" + "sample_" + b + ".2.fq" , "a")
		
	total_rec = []
	drop_no_bc = []
	drop_low_q = []
		
	for k in sorted(barcodes_dict.keys(), key = int)[::-1]: # look for the short ones first	
	#	with gzip.open(outfile_firstreads, "w") as OUTF:
	#		with gzip.open(first_reads, "r") as INF:
		discard1 = open(first_reads + ".discard", "a")
		discard2 = open(second_reads + ".discard", "a")
		print "checking barcodes of length", k
		logfile.write("# barcodes of length " + str( k) + "\n")
		total_entry_count = 0
		successful_demultiplexed = 0
		no_barcode = 0
		low_q = 0
		with open(first_reads + ".tmp", "r") as INF1:
			with open(second_reads + ".tmp", "r") as INF2:
				for [name1, seq1, qual1], [name2, seq2, qual2] in zip(readfq(INF1),readfq(INF2)):
					total_entry_count += 1
					qual1_numeric = [ ord(x) - 33 for x in qual1.rstrip('\n') ] ## from http://pythonforbiologists.com/index.php/business-card-fastq-parser-exercise/
					qual2_numeric = [ ord(x) - 33 for x in qual2.rstrip('\n') ]
					if len( [ x for x in qual1_numeric[:8] if x >= 20] ) == 8:
						good = True					
						i = 0
						while i < readlen-window_size-1:
							i += 1
							if not sum(qual1_numeric[i:i+window_size+1]) >= min_qual_sum:
								good = False
								break
							if not sum(qual2_numeric[i:i+window_size+1]) >= min_qual_sum:
								good = False
								break
						if not seq2.startswith(args.rmot2):
							good = False
						if good == True:
							bc_seen = None
							for b in barcodes_dict[k]:
								if seq1.startswith(b + args.rmot1):					
									outentry1 = "@" + name1 + "/1" + "\n" + seq1[k:] + "\n+\n" + qual1[k:] + "\n"
									outentry2 = "@" + name2 + "/2" + "\n" + seq2 + "\n+\n" + qual2 + "\n"
									output_files_dict[b+"1"].write( outentry1 )
									output_files_dict[b+"2"].write( outentry2 )
									successful_demultiplexed += 1
									bc_seen = True
									break
							if not bc_seen:
								outentry1 = "@" + name1 + "\n" + seq1 + "\n+\n" + qual1 + "\n"
								outentry2 = "@" + name2 + "\n" + seq2 + "\n+\n" + qual2 + "\n"
								discard1.write( outentry1 )
								discard2.write( outentry2 )	
								no_barcode += 1							
						else:
#							discard1.write( outentry1 )
#							discard2.write( outentry2 )
							low_q += 1
					else:
#						discard1.write( outentry1 )
#						discard2.write( outentry2 )
						low_q += 1
						
		discard1.close()
		discard2.close()
		os.rename(first_reads + ".discard", first_reads + ".tmp")
		os.rename(second_reads + ".discard", second_reads + ".tmp")			
		logfile.write( "total read pairs:		" + str(total_entry_count) + "\n" + "passed:		" + str(successful_demultiplexed) + "\n" + "dropped:	" + str(low_q + no_barcode) + "\n" + "prop dropped:	" + str( float(low_q + no_barcode)/float( total_entry_count ) ) + "\n" +
		"no barcode:	" + str(no_barcode) + "\n" )
		total_rec.append( total_entry_count )
		drop_no_bc.append( no_barcode )
		drop_low_q.append( low_q )
		
	os.rename(first_reads + ".tmp", first_reads + ".final_discards")			
	os.rename(second_reads + ".tmp", second_reads + ".final_discards")
	tot_read_pairs = max(total_rec)
	final_discards = drop_no_bc[-1] + drop_low_q[0]

	logfile.write( "########\n" +
	"overall total:		" + str(tot_read_pairs) + "\n" +
	 "overall passed:		" + str(tot_read_pairs - final_discards ) + "\n" +
	 "prop. overall passed:		" + str( float(tot_read_pairs - final_discards)/ tot_read_pairs ) + "\n" +
	 "no barcode:		" + str( drop_no_bc[-1] ) + "\n" +
	 "low quality:		" + str( drop_low_q[0] ) + "\n" )
	
	print "Done!"
	logfile.close()
	