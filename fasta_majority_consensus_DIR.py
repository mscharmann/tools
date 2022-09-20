
import sys, os


def read_fasta (infile):

	fasta_dict = {}
	seq = ""
	name = "dummy"
	with open(infile, "r") as F:
		for line in F:
			if line.startswith(">"):
				fasta_dict[name] = seq
				name = line.lstrip(">").rstrip("\n")
				seq = ""
			else:
				seq += line.strip("\n")

	# last record:
	fasta_dict[name] = seq
	del fasta_dict["dummy"]
	seqs = list( fasta_dict.values() )

	return seqs



def make_maj_consensus (seqs):
	
	consensus_seq = ""
	for i in range(len(seqs[0])):
		NUCs = [x[i] for x in seqs if not x[i] == "-"]
		unique_NUCs = list(set(NUCs))
		counts = []
		for A in unique_NUCs:
			counts.append( NUCs.count(A) )
		maxidx = counts.index(max(counts))
		consensus_seq += unique_NUCs[maxidx]
	return consensus_seq 
	
####

if len(sys.argv) != 2:
	print ( "usage: python translate_CDS_to_PEP.py cds_fasta_dir" )
	exit()

infiles = [x for x in os.listdir(sys.argv[1]) if x. endswith(".CDS.aln")]

outlines = []
cnt = 0
for f in infiles:
	cnt += 1
	print(cnt)
	seqs = read_fasta(sys.argv[1] + "/" + f)
	if len(seqs) > 0:
		cons = make_maj_consensus(seqs)
		outlines.append(">" + f.strip(".CDS.aln"))
		outlines.append(cons)
	else:
		with open("maj_cons_problems_log.txt", "a") as O:
			O.write(f + "\n")
			print ( f )
	
with open("maj_cons.fasta", "w") as O:
	O.write("\n".join(outlines) + "\n")

	
	


