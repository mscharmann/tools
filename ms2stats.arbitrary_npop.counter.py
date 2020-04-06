#!/usr/local/bin/python
# Python 2.7
# 28 Feb 2018
#
# ms2stats.arbitrary_npop.counter.py
#
# Mathias Scharmann

# INPUTS: 	- stdin in ms format
#			- "spinput.txt"-file: specifies number of loci, sample sizes and site-lengths for each locus

# OUTPUT: 	- "ABCstat.txt"-file with selected popgen summary statistics

# example usage:
"""
just counts ms output and sends to stdout

python $HOME/tools/draw_ms_priors.1.py -bpfile $bpfile -argfile $argfile  | $toolbase/msnsam/msnsam tbs $msnsam_arg2 $msnsam_model_cmd | python ms2stats.counter.py | python ms2stats.stats.py


"""

##############
# HEAD
##############

import sys

##############
# FUNCTIONS
##############

def determine_npop_from_spinput ():
	
	line_cnt = 0
	
	with open("spinput.txt", "r") as INFILE:
		INFILE.readline() # discard first line (empty)
		n_loci_per_replicate = int( INFILE.readline().strip("\n") )
		for line in INFILE:
			if len(line) > 1:
				# count non-empty lines
				line_cnt += 1
	npop = ( float(line_cnt-2)/n_loci_per_replicate ) -1
	if not npop.is_integer():
		print "failed to determine npop from spinput.txt"
		exit()
	return int(npop)

def read_spinput(npop):
	
	nsam_locilengths = []
	
	with open("spinput.txt", "r") as INFILE:
		INFILE.readline() # discard first line (empty)
		n_loci_per_replicate = int( INFILE.readline().strip("\n") )
		
		loci_cnt = 0
		while loci_cnt < n_loci_per_replicate:
			loc_sampling_scheme = []
			for p in range(npop):				
				nsam_p = int( INFILE.readline().strip("\n") )
				loc_sampling_scheme.append( nsam_p )
			locus_length = int( INFILE.readline().strip("\n") )
			loc_sampling_scheme.append( locus_length )
			nsam_locilengths.append( loc_sampling_scheme )
			
			loci_cnt += 1
		
		n_replicates = int( INFILE.readline().strip("\n") )
	
	return n_replicates, nsam_locilengths
	
def ms2counts(ms_out, nsam_and_locuslength):


	# catch loci without segsites but return fixed counts!
	try:

		sites_cnts = []
		nsites = len( ms_out[0] )

		for s in range( nsites ):
			totalsam = 0
			s_counts = []
			for pop in range(len(nsam_and_locuslength)-1):
				nsam = nsam_and_locuslength[pop]
	#			print "nsam", nsam 

				pop_samples = ms_out[totalsam:totalsam+nsam]
		#		print pop_samples
				totalsam += nsam
			
				resh = [ [row[idx] for row in pop_samples]  for idx in range(len(pop_samples[0])) ] 
	#			print len(resh)
	
				cnt1 = resh[s].count("1")
				cnt0 = nsam - cnt1
				s_counts.append( cnt1 )
				s_counts.append( cnt0 )
			sites_cnts.append( s_counts )
	
	except IndexError:

		fixed = [ ]
		for pop in range(len(nsam_and_locuslength)-1):
			nsam = nsam_and_locuslength[pop]
			fixed.append(nsam)
			fixed.append(0)
		sites_cnts = [ fixed ]
		
#	print "newmethod", sites_cnts
	"""
	final structure of output is:
	site_1_pop_1_count_1 site_1_pop_1_count_0 site_1_pop_2_count_1 site_1_pop_2_count_0 site_1_pop_3_count_1 site_1_pop_3_count_0
	site_2_pop_1_count_1 site_2_pop_1_count_0 site_2_pop_2_count_1 site_2_pop_2_count_0 site_2_pop_3_count_1 site_2_pop_3_count_0
	site_3_pop_1_count_1 site_3_pop_1_count_0 site_3_pop_2_count_1 site_3_pop_2_count_0 site_3_pop_3_count_1 site_3_pop_3_count_0
	"""
	return sites_cnts


##############
# MAIN
##############

if __name__ == "__main__":
	
	npop = determine_npop_from_spinput ()

	n_replicates, nsam_locilengths = read_spinput(npop) 

	n_loci = len( nsam_locilengths )

	locilengths = [ x[-1] for x in nsam_locilengths ]

	locus_cnt = 0
	replicate_cnt = 0
	
	slinecnt = 0
	
	samples = []
	loci_counts_of_replicate = []
	
	preceding_line = "placeholderstring"
	
	outlines = []

			
	# discard first 3 lines:
 	sys.stdin.readline()
 	sys.stdin.readline()
 	sys.stdin.readline()
# 	first_locus_cmd = sys.stdin.readline().split("	")
# 	# nsam_locuslength = [ int(first_locus_cmd[idx]) for idx in [5,6,4] ]
# 	nsam_locuslength = nsam_locilengths[ 0 ]
# 	locilengths.append( nsam_locuslength[2] )
		
	for line in sys.stdin:
#		print line
		
		try:
			first_char = int( line[0] ) # this means it is a sample line (most lines are)
			samples.append( line.strip("\n") )
			slinecnt += 1
#			print "read chrom", slinecnt
    		
		except ValueError: # first char is not an int, so a character line
			
			if line.startswith("//"):
				# this means that a new locus starts
				locus_cnt += 1
				samples = []
				# nsam_locuslength = [ int( line.split("	")[idx] ) for idx in [5,6,4] ]
				# locilengths.append( nsam_locuslength[2] )	
									
			elif len(line) == 1:	
				# empty line : means the previous locus is finished OR that there were no samples in this locus because of segsites: 0
				if len(preceding_line) == 1:
					# skip, because locus was already finished with preceding empty line
					continue				
					
				else:
#					print "locus_finished, counting", len(samples)
			
					nsam_locuslength = nsam_locilengths[locus_cnt-1] # zero-based index!
					site_cnts = ms2counts(samples, nsam_locuslength)
					
					# and send to stdout
					sys.stdout.write( " ".join( [ " ".join( [str(cnt) for cnt in SNP] ) for SNP in site_cnts] ) + "\n")

					if (locus_cnt ) == n_loci :
						# this means that one replicate == set of x loci is full and we have to harvest it
# 						print "replicate_finished, harvesting stats"
# 						print locus_cnt, len(loci_counts_of_replicate)
						
						sys.stdout.write( "\n" )

						locus_cnt = 0
						
			else:
				continue
		preceding_line = line			
	
	# end of stdin: if stdin does not end with an empty line (which is the regular ms output), replicate harvest is not triggered! So do it here.
	if len(line) > 1:
		
		nsam_locuslength = nsam_locilengths[locus_cnt-1] # zero-based index!
		site_cnts = ms2counts(samples, nsam_locuslength )
		sys.stdout.write( " ".join( [ " ".join( [str(cnt) for cnt in SNP] ) for SNP in site_cnts] ) + "\n")
		
	# always end stdout with an empty line that triggers harvest of last replicate for the stats module
	sys.stdout.write( "\n")
			
