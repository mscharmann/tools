
import sys

for line in sys.stdin:
	outl = []
	for field in line.strip("\n").split("\t"):
		if field.startswith("./."):
			outl.append(".")
		elif field.startswith(".:."):
			outl.append(".")
		else:
			outl.append( field )
	sys.stdout.write( "\t".join( outl ) + "\n" )	
		

