import sys
import os


# - each pairwise alignment must have at least 10 codons = 30 bases
# - each pairwise alignment must have at least 3 SNPs
# - eventual stop codons are dropped.
# if conditions not met, omega is returned as NA



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
	return fasta_dict



def export_clean_fasta (seqdatadict):
	
	
	stop_codons = set(["TAG","TGA","TAA"])
	
	try:
		[seq1,seq2] = list(seqdatadict.values())
	except ValueError:
		# there was only one sequence in the fasta; this can happen if prequal drops one seq entirely...
		return None
		
	outseq1 = ""
	outseq2 = ""
	
	# first get rid of stop codons
	for triplet in range(0, len(seq1), 3):
		if not seq1[triplet:triplet+3] in stop_codons:
			if not seq2[triplet:triplet+3] in stop_codons:
				outseq1 += seq1[triplet:triplet+3]
				outseq2 += seq2[triplet:triplet+3]
				
	seq1 = outseq1
	seq2 = outseq2	
	outseq1 = ""
	outseq2 = ""

	snps = 0
	for pos in range(len(seq1)):
		if not seq1[pos] == "-":
			if not seq2[pos] == "-":
				outseq1 += seq1[pos]
				outseq2 += seq2[pos]
				if seq1[pos] != seq2[pos]:
					snps += 1
		
	if len(outseq1) > 30 and snps >= 3:	
		outlines = []
		outlines.append(">sequence_1")
		outlines.append(outseq1)
		outlines.append(">sequence_2")
		outlines.append(outseq2)
		outfilename = "clean_for_yn00.fasta"
		with open(outfilename, "w") as F:
			F.write("\n".join(outlines)+"\n")
		return [outfilename, len(outseq1), snps]		
	else:
		return [None, "<=30", "<3"]


def run_PAML_codeml_yn00 (infile):
	rundirname = "rundir_" + infile
	try:
		os.mkdir(rundirname)
	except OSError:
		None
	ctlstring = """
	  seqfile = ../{} * sequence data file name
      outfile = yn.txt           * main result file
      verbose = 0  * 1: detailed output (list sequences), 0: concise output

        icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below

    weighting = 0  * weighting pathways between codons (0/1)?
   commonf3x4 = 0  * use one set of codon freqs for all pairs (0/1)? 
*       ndata = 1	
	""".format(infile)
	with open(rundirname + "/yn00.ctl", "w") as O:
		O.write(ctlstring)
	cmd = """
	cd {}
	yn00 yn00.ctl &>/dev/null	
	""".format(rundirname)
	os.system(cmd)
	
	#omega = "NA"
	
	## parse the output file:
	with open(rundirname + "/yn.txt", "r") as INF:
		for line in INF:
			fields = line.strip("\n").split()
			if len(fields) > 0:
				if fields[0] == "2":
					omega = fields[6]
					dN = fields[7]
					dN_SE = fields[9]
					dS = fields[10]
					dS_SE = fields[12]
					break
	os.system("rm -r " + rundirname)
	os.system("rm " + infile)
	
	return [omega, dN, dN_SE, dS, dS_SE]
	

############## MAIN

#treefile = sys.argv[1]

all_infiles = [x for x in os.listdir(sys.argv[1]) if x.endswith(".CDS.aln")]
print(all_infiles)


OUTF = open("yn00.report.txt", "w")

OUTF.write("# PAML yn00\n")
OUTF.write("# method: Yang Z, Nielsen R (2000) Estimating synonymous and nonsynonymous substitution rates under realistic evolutionary models. Mol. Biol. Evol. 17:32-43\n")
OUTF.write("# alignment" + "\t" + "\t".join(["aln_length","substitutions","omega", "dN", "dN_SE", "dS", "dS_SE"]) + "\n")

for aln in all_infiles:
	print (aln)
	seqdatadict = read_fasta (sys.argv[1] + "/" + aln)
	[fasta_name, aln_len, snps] = export_clean_fasta (seqdatadict)
	if fasta_name:
		[omega, dN, dN_SE, dS, dS_SE] = run_PAML_codeml_yn00 (fasta_name)
	else:
		[omega, dN, dN_SE, dS, dS_SE] = ["NA"]*5
	OUTF.write(aln + "\t" + "\t".join([str(aln_len),str(snps),omega, dN, dN_SE, dS, dS_SE]) + "\n")

OUTF.close()

		
		
		
