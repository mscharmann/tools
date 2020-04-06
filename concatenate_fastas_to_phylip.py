"""
Filter the aligned and cleaned ortholog matrices by number of taxa and characters
Write out name of matrices that passed the filter
Also write out supermatrix stats
MS mod: sleep a few seconds after phyutility writes its outout, to make sure the file is readable in the next step
"""
from seq import read_fasta_file
from tree_utils import get_name
import sys,os,time

def concatenate(clnDIR,numofsitesFilter,numoftaxaFilter,seqtype,outname):
	"""filter cleaned alignments and concatenate"""
	if clnDIR[-1] != "/": clnDIR += "/"
	sites_filter = int(numofsitesFilter)
	taxa_filter = int(numoftaxaFilter)
	assert seqtype == "aa" or seqtype == "dna","Input data type: dna or aa"
	model = "WAG"  if seqtype == "aa" else "DNA"
	
	print "Filtering ortholog matrixes"
	selected = [] # list of alignment file names that pass the filters
	for i in os.listdir(clnDIR):
		if i.endswith(".aln-cln"):
			seqlist = read_fasta_file(clnDIR+i)
			num_seq = len(seqlist)
			num_col = len(seqlist[0].seq)
			if num_seq >= taxa_filter and num_col >= sites_filter:
				selected.append(i)
	print len(selected),"matrices passed the filter"

	print "Getting matrix occupancy stats"
	taxon_occupancy = {}
	#key is taxon name, value is [times present in a matrix,total length for this taxon]
	total_aligned_len = 0 #record how long the final concatenated matrix is
	if seqtype == "aa":
		cmd = "phyutility -concat -aa -out "+outname+".nex -in "
	else: cmd = "phyutility -concat -out "+outname+".nex -in "
	for i in selected:
		cmd += clnDIR+i+" "
		seqlist = read_fasta_file(clnDIR+i)
		total_aligned_len += len(seqlist[0].seq)
		for s in seqlist:
			taxonid = get_name(s.name)
			if taxonid not in taxon_occupancy:
				taxon_occupancy[taxonid] = [0,0]
			taxon_occupancy[taxonid][0] += 1
			taxon_occupancy[taxonid][1] += len((s.seq.replace("-","")).replace("?",""))
	cmd += "\n"
	
	total_ortho = len(selected)
	with open(outname+"_taxon_occupancy_stats","w") as outfile:
		outfile.write("taxon\t#orthologs\t#total_charactors\tperc_orthologs\tperc_charactors\n")
		sum_char = 0
		for taxon in taxon_occupancy:
			times,chars = taxon_occupancy[taxon][0],taxon_occupancy[taxon][1]
			sum_char += chars
			out = taxon+"\t"+str(times)+"\t"+str(chars)+"\t"
			out += str(times/float(total_ortho))+"\t"+str(chars/float(total_aligned_len))+"\n"
			outfile.write(out)
		total_taxa = len(taxon_occupancy)
		out = "\nSupermatrix dimension "+str(total_taxa)+" taxa, "
		out += str(total_ortho)+" loci and "+str(total_aligned_len)+" aligned columns\n"
		out += "Overall matrix occupancy "+str(sum_char/float(total_taxa*total_aligned_len))+"\n"
		outfile.write(out)
	
	print "Supermatrix taxon occupancy stats written to",outname+"_taxon_occupancy_stats"
	print "Waiting for concatenation to finish. This may take several minutes..."
	
	phy_dict = {}
	for tax in taxon_occupancy.keys():
		phy_dict[tax] = []

	left_idx = 1
	right_idx = 0
	modelfile_lines = []
	for i in selected:
		mline = model + "," + i + "=" + str( left_idx ) + "-" + str( right_idx )		
		seqlist = read_fasta_file(clnDIR+i)
		aligned_len = len(seqlist[0].seq)
		right_idx += aligned_len		
		mline = model + "," + i + "=" + str( left_idx ) + "-" + str( right_idx )	
		left_idx += aligned_len
		print mline
		modelfile_lines.append( mline )
		fasta_dict = {s.name : s.seq for s in seqlist }
		for tax in taxon_occupancy.keys():
			try:
				phy_dict[tax].append( fasta_dict[tax] )
			except KeyError:
				phy_dict[tax].append("-"*aligned_len)
	
#	print phy_dict 

	with open(outname+".model","w") as outfile2:
		outfile2.write( "\n".join( modelfile_lines ) )	
	
	with open(outname + ".phy", "w") as OUTFILE:
		ntaxa = len(phy_dict.keys())
		len_align = len( "".join( phy_dict[phy_dict.keys()[0]] ) ) # gets value of first element in dictionary -> this OK since all seqs have same length
		header = str(ntaxa) + " " + str(len_align) + "\n"
		OUTFILE.write(header)
		for sample, seq in phy_dict.items():
			out_line = sample + "    " + "".join(seq) + "\n"
			OUTFILE.write(out_line)
	OUTFILE.close()

	print "outfiles written",outname+".phy",outname+".model"


if __name__ == "__main__":
	if len(sys.argv) != 6:
		print "usage: python concatenate_matrices.py aln-clnDIR numofsitesFilter numoftaxaFilter dna/aa outname"
		sys.exit()
		
	clnDIR,numofsitesFilter,numoftaxaFilter,seqtype,outname = sys.argv[1:]
	concatenate(clnDIR,numofsitesFilter,numoftaxaFilter,seqtype,outname)
