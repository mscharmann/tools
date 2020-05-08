#!/usr/local/bin/python
# Python 2.7.6 / 2.6
# phylip2stockholm.py
# 20 Feb 2015
# Mathias Scharmann

"""

usage:

toolbox=/gdc_home3/schamath/tools
python $toolbox/phylip2stockholm.py raff_complex_Brunei_r0.75.phy

"""

# prob. best in python 2.6 since Bio is nor installed for 2.7


from Bio import AlignIO
import sys


infile=sys.argv[1]
outfile=infile.rstrip(".phy")+".sth"


input_handle = open(infile, "rU")
output_handle = open(outfile, "w")
 
alignments = AlignIO.parse(input_handle, "phylip-relaxed")
AlignIO.write(alignments, output_handle, "stockholm")
 
output_handle.close()
input_handle.close()

print "made stockholm from phylip,\nDone!"
