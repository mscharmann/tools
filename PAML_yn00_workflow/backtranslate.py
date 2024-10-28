
import sys

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


def cds_aln_from_aa_aln_and_cds_unalign ( aa_aln_file, cds_unalign_infile ):
	
	aa_seqs = read_fasta( aa_aln_file )
	cds_unaligned_seqs = read_fasta( cds_unalign_infile )
	
	aln_vs_raw_dict = {}
	for id, seq in aa_seqs.items():
		aln_vs_raw_dict[ id ] = [ seq ]
	for id, seq in cds_unaligned_seqs.items():
		try:
			aln_vs_raw_dict[ id ].append( seq )
		except KeyError:
			None
	
	out_lines = []
	for name, v in aln_vs_raw_dict.items():
		aa_aln = v[0]
		cds_unaln = v[1]
		cds_aln = ""
		cds_idx = 0
		for idx, char in enumerate( aa_aln ):
			if char != "-":
				cds = cds_unaln[ cds_idx : cds_idx+3 ]
				cds_idx += 3
			else:
				cds = "---"
			# print idx, char, cds
			cds_aln += cds
		out_lines.append( ">" + name )
		out_lines.append( cds_aln )
	
	outname = aa_aln_file + ".backtranslated.aln"
	with open( outname, "w" ) as OUTF:
		OUTF.write(	"\n".join( out_lines ) + "\n")
	
####

if len(sys.argv) != 3:
	print ("usage: python backtranslate.py aa_aln_file cds_unalign_infile")
	exit()
	
cds_aln_from_aa_aln_and_cds_unalign(sys.argv[1],sys.argv[2])
