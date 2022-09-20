##
## merges two tab-separated files based on the strings in a column, fills NA if line not exists!

import sys

infile1 = sys.argv[1]
column_1 = sys.argv[2]
infile2 = sys.argv[3]
column_2 = sys.argv[4]



with open(infile1, "r") as INF:
	indict1 = { line.strip("\n").split("\t")[column_1] : line.strip("\n") for line in INF }

with open(infile2, "r") as INF:
	indict2 = { line.strip("\n").split("\t")[column_2] : line.strip("\n") for line in INF }

all_keys = set( indict1.keys() + indict2.keys() )

outlines = []
for k in all_keys:	
	outl = ""	
	try:
		outl += indict1[k]
	except KeyError:
		outl += "\t".join( ["NA"]*len( indict1[indict1.keys()[0]] ) )
	outl += "\t"	
	try:
		outl += indict2[k]
	except KeyError:
		outl += "\t".join( ["NA"]*len( indict2[indict2.keys()[0]] ) )
	outlines.append( outl )
	

sys.stdout.write("\n".join(outlines))

