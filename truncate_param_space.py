

import sys

tau_threshold = 0.025

#######################################
## read
with open("parameters.mig0_t5_p0.flexighosts.txt", "r") as INF:
	params = [x for x in INF]


with open("ABCstat.mig0_t5_p0.flexighosts.txt", "r") as INF:
	stats = [x for x in INF]


# WE NEED TO REMOVE ANY TIMEs < 0.0025

outlines_params = []
outlines_stats = []
for param_line, stats_line in zip( params, stats ):
	try:	
		if float( param_line.split("\t")[ 6 ] ) >= tau_threshold and float( param_line.split("\t")[ 7 ] ) >= tau_threshold:
			outlines_params.append( param_line )
			outlines_stats.append( stats_line )
	except ValueError:
		outlines_params.append( param_line )
		outlines_stats.append( stats_line )
		
print "retained:	", len(outlines_params)-1			
		
with open("parameters.mig0_t5_p0.flexighosts.txt_filtered", "w") as OUTF:
	OUTF.write("".join(outlines_params))

with open("ABCstat.mig0_t5_p0.flexighosts.txt_filtered", "w") as OUTF:
	OUTF.write("".join(outlines_stats))


#######################################
## read
with open("parameters.mig0_t5_p0.txt", "r") as INF:
	params = [x for x in INF]


with open("ABCstat.mig0_t5_p0.txt", "r") as INF:
	stats = [x for x in INF]


outlines_params = []
outlines_stats = []
for param_line, stats_line in zip( params, stats ):
	try:	
		if float( param_line.split("\t")[ 5 ] ) >= tau_threshold and float( param_line.split("\t")[ 5 ] ) >= tau_threshold:
			outlines_params.append( param_line )
			outlines_stats.append( stats_line )
	except ValueError:
		outlines_params.append( param_line )
		outlines_stats.append( stats_line )
		
print "retained:	", len(outlines_params)-1			
		
with open("parameters.mig0_t5_p0.txt_filtered", "w") as OUTF:
	OUTF.write("".join(outlines_params))

with open("ABCstat.mig0_t5_p0.txt_filtered", "w") as OUTF:
	OUTF.write("".join(outlines_stats))



#######################################
## read
with open("parameters.mig20_t5_p0.95.txt", "r") as INF:
	params = [x for x in INF]


with open("ABCstat.mig20_t5_p0.95.txt", "r") as INF:
	stats = [x for x in INF]


outlines_params = []
outlines_stats = []
for param_line, stats_line in zip( params, stats ):
	try:	
		if float( param_line.split("\t")[ 5 ] ) >= tau_threshold and float( param_line.split("\t")[ 5 ] ) >= tau_threshold:
			outlines_params.append( param_line )
			outlines_stats.append( stats_line )
	except ValueError:
		outlines_params.append( param_line )
		outlines_stats.append( stats_line )


print "retained:	", len(outlines_params)-1		
		
		
with open("parameters.mig20_t5_p0.95.txt_filtered", "w") as OUTF:
	OUTF.write("".join(outlines_params))

with open("ABCstat.mig20_t5_p0.95.txt_filtered", "w") as OUTF:
	OUTF.write("".join(outlines_stats))

