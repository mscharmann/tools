import sys


if len(sys.argv) != 3:
	print "usage: python make_arbitrary_model_file.py phylip_file length_of_partitions (integer)"
	exit()



	
pyhlip_in = sys.argv[1]
plength = int( sys.argv[2] )


with open(pyhlip_in, "r") as INF:
	ncol = int( INF.readline().strip("\n").split()[1] )



print ncol
startcol = 1
endcol = 0
cnt = 0

partitions = []
while endcol < ncol:
	cnt += 1
	endcol += plength
	outl = "DNA,partition_{0}={1}-{2}".format(cnt, startcol, endcol  )
	partitions.append(outl)
	startcol += plength
	

cnt += 1
endcol = ncol
outl = "DNA,partition_{0}={1}-{2}".format(cnt, startcol, endcol  )
partitions.append(outl)


print "made {0} partitions of length {1}".format( len(partitions), plength )

with open(pyhlip_in + ".model.txt", "w") as OUTF:
	OUTF.write( "\n".join( partitions ) )