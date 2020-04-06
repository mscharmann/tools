import sys

#Usage: filter_fasta_per_ids.py input.fasta filter_ids.txt output.fasta

input_file =sys.argv[1]
id_file =sys.argv[2]
output_file =sys.argv[3]
wanted = set(line.rstrip("\n").split(None,1)[0] for line in open(id_file))
print("Found %i unique identifiers in %s" % (len(wanted), id_file))

input_fasta = {}
with open(input_file, "r") as INFILE:
	for line in INFILE:
		if line.startswith(">"):
			id = line.strip(">").strip("\n")
		elif line.startswith("N"):
			input_fasta[id] = line.strip("\n")

coreids = {}
for id in wanted:
	coreids[ id.split("_")[0] ] = id

the_order = sorted(coreids.keys(), key = int)	
the_order = [ coreids[x] for x in the_order[:] ]

print the_order
			
outlines = []
for id in the_order:
		outlines.append( ">"+id+"\n"+input_fasta[id] )

with open(output_file, "w") as OUTPUT:
	OUTPUT.write("\n".join(outlines))
