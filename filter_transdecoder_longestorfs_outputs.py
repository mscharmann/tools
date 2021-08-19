#### filter_transdecoder_longestorfs_outputs.py

"""
python filter_transdecoder_longestorfs_outputs.py DIR

- created a new directory DIR.cleaned with the outputs, ready for TransDecorder Predict
- retains any ORFs >= 100 AA
- ORFs < 100 AA are retained IF they have a BLAST hit

"""

import os, sys

DIR = sys.argv[1]
if not DIR.endswith("/"):
	DIR = DIR + "/"

outdir = DIR.rstrip("/") + ".cleaned"
try:
	os.mkdir(outdir)	
except OSError:
	print "Warning: output DIR already exists, will overwrite"

pepfile = "longest_orfs.pep"
cdsfile = "longest_orfs.cds" 
blastfile = "blastp.outfmt6"
gff_file = "longest_orfs.gff3"

def read_fasta (infile):
	
	seq = ""
	id = "dummy"
	out_dict = {}
	with open(infile, "r") as F:
		for line in F:
			if len(line) > 1:
				if line.startswith(">"):
					out_dict[id] = seq
					seq = ""
					id = line.strip("\n")
				else:
					seq += line.strip("\n")
	out_dict[id] = seq
	del out_dict["dummy"]
	return out_dict


print "preparing"
blast_lines = {}
with open(DIR + "/" + blastfile, "r") as F:
	for line in F:
		if len(line) > 1:
			fields = line.split("\t")
			blast_lines[fields[0]] = line


seqs_to_retain = []
all_blast_hit_ids = set( blast_lines.keys() )
with open(DIR + pepfile, "r") as PEPFILE:
	for line in PEPFILE:
		if line.startswith(">"):
			id = line.lstrip(">").split(" ")[0]
			AA_length = int( line.split("len:")[1].split(" ")[0] )
			if AA_length >= 100:
				seqs_to_retain.append(id)
			else:
				if id in all_blast_hit_ids:
					seqs_to_retain.append(id)
			

print len( seqs_to_retain )
seqs_to_retain = set(seqs_to_retain)
print len( seqs_to_retain )



print "parsing PEP file"
pepseqs = read_fasta (DIR + pepfile)
outl = [x + "\n" + pepseqs[x] + "\n" for x in pepseqs.keys() if x.strip(">").split(" ")[0] in seqs_to_retain]
with open(outdir + "/" + pepfile, "w") as OUTF:
	OUTF.write("".join(outl))


print "parsing CDS file"
cdsseqs = read_fasta (DIR + cdsfile)
outl = [x + "\n" + cdsseqs[x] + "\n" for x in cdsseqs.keys() if x.strip(">").split(" ")[0] in seqs_to_retain]
with open(outdir + "/" + cdsfile, "w") as OUTF:
	OUTF.write("".join(outl))


print "parsing blast file"
outl = [ ]
for x in seqs_to_retain:
	try:
		outl.append( blast_lines[x] )
	except KeyError:
		None

with open(outdir + "/" + blastfile, "w") as OUTF:
	OUTF.write("".join(outl))


# features = ['five_prime_UTR', 'exon', 'mRNA', 'CDS', 'three_prime_UTR', 'gene']
print "parsing gff"
outl = []
with open(DIR + gff_file, "r") as F:
	for line in F:
		if not line.startswith("#"):
			if len(line) > 1:
				fields = line.split("\t")
				if fields[2] == 'gene' or fields[2] == 'mRNA':
					id = line.split("~~")[1].split(";")[0]
				else:
					id = line.strip("\n").split("Parent=")[1].split(";")[0]
				if id in seqs_to_retain:
					outl.append(line)
		else:
			outl.append(line) 


with open(outdir + "/" + gff_file, "w") as OUTF:
	OUTF.write("".join(outl))


print "copying base_freqs.dat"
with open(DIR + "base_freqs.dat", "r") as F:
	with open(outdir + "/" + "base_freqs.dat", "w") as O:
		O.write(F.read())


#####################################

