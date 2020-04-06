
import sys

"""
takes a phylip alignment of codons (!) and extracts 4-fold degenerate sites, returns as phylip format!

the aim:

Estimation of Evolutionary Distance.

as in:
Kumar, S. and Subramanian, S. (2002). Mutation rates in mammalian genomes. PNAS 99, 803-808.

methods section:
Evolutionary divergence (d4) between sequences at fourfold-degenerate sites was estimated by using the Tamura-Nei method (22) to correct for multiple hits by accounting for transition/transversion rate and base-frequency biases. 

this "standardised" measure of evolutionary distance will allow to place the Nepenthes radiation into comparison with familiar evolutionary splits (e.g. human - chimp) and the divergence in other radiations (e.g. mammals)

## but maybe not necessary? Here from chimp genome paper:

The Chimpanzee Sequencing and Analysis Consortium (2005). Initial sequence of the chimpanzee genome and comparison with the human genome. Nature 437, 69-87.

comparing ALL sites / whole genomes, coding or not:
"Single-nucleotide substitutions occur at a mean rate of 1.23% between copies of the human and chimpanzee genome, with 1.06% or less corresponding to fixed divergence between the species."


## results from ortho121.noDIONtrim.supermatrix.8.phy ::
EPHI$$$RAFF [118015, 7759914, 1.52]
PERV$$$RAFF [295282, 7275420, 4.06]
GRAC$$$JACQ [123896, 7444206, 1.66]
HAMA$$$VOGE [101816, 7679949, 1.33]
PERV$$$RAJA [297379, 7238562, 4.11]
DION$$$HEMS [738745, 3584850, 20.61]
HEMS$$$KHAS [177653, 7698222, 2.31]
GRAC$$$KHAS [183544, 7530744, 2.44]
HEMS$$$VOGE [115971, 7815444, 1.48]
KHAS$$$VOGE [184881, 7571886, 2.44]
AMPU$$$DION [673302, 3278175, 20.54]
EPHI$$$RAJA [71260, 7741446, 0.92]
AMPU$$$SPGI [85673, 7376475, 1.16]
KHAS$$$SPGI [180948, 7659516, 2.36]
DION$$$KHAS [732377, 3557079, 20.59]
DION$$$RAJA [740594, 3581154, 20.68]
EPHI$$$GRAC [124757, 7605990, 1.64]
EPHI$$$HEMS [116531, 7803177, 1.49]
DION$$$EPHI [725690, 3507096, 20.69]
JACQ$$$RAFF [116941, 7591113, 1.54]
AMPU$$$RAJA [110000, 7341585, 1.5]
HAMA$$$SPGI [123300, 7765518, 1.59]
HEMS$$$SPGI [33608, 7962600, 0.42]
HAMA$$$KHAS [188053, 7572234, 2.48]
DION$$$VOGE [735265, 3559587, 20.66]
SPGI$$$VOGE [119838, 7775157, 1.54]
RAFF$$$SPGI [34808, 7920009, 0.44]
AMPU$$$JACQ [112603, 7118274, 1.58]
JACQ$$$RAJA [84881, 7574754, 1.12]
AMPU$$$EPHI [113507, 7254504, 1.56]
RAJA$$$VOGE [70434, 7772220, 0.91]
HAMA$$$HEMS [120623, 7826595, 1.54]
EPHI$$$SPGI [120046, 7747797, 1.55]
HAMA$$$JACQ [104090, 7519470, 1.38]
DION$$$RAFF [738448, 3583872, 20.6]
JACQ$$$KHAS [183314, 7407702, 2.47]
GRAC$$$PERV [297404, 7150977, 4.16]
AMPU$$$VOGE [112877, 7250490, 1.56]
EPHI$$$KHAS [184819, 7518780, 2.46]
RAFF$$$VOGE [117296, 7777722, 1.51]
GRAC$$$HAMA [127777, 7638201, 1.67]
AMPU$$$KHAS [170073, 7150320, 2.38]
AMPU$$$HAMA [116948, 7278111, 1.61]
HAMA$$$RAJA [98273, 7755666, 1.27]
HEMS$$$JACQ [115613, 7623225, 1.52]
EPHI$$$PERV [297881, 7140855, 4.17]
DION$$$HAMA [734129, 3545382, 20.71]
JACQ$$$VOGE [87684, 7517106, 1.17]
JACQ$$$SPGI [119539, 7585146, 1.58]
PERV$$$VOGE [297955, 7177452, 4.15]
EPHI$$$JACQ [87773, 7484832, 1.17]
GRAC$$$VOGE [124401, 7648770, 1.63]
HEMS$$$RAJA [112612, 7909185, 1.42]
RAFF$$$RAJA [114180, 7868004, 1.45]
GRAC$$$SPGI [108833, 7721205, 1.41]
DION$$$PERV [714025, 3465399, 20.6]
HAMA$$$PERV [301968, 7196019, 4.2]
KHAS$$$PERV [283105, 7161201, 3.95]
HEMS$$$RAFF [28144, 8006055, 0.35]
DION$$$SPGI [740762, 3589836, 20.63]
PERV$$$SPGI [296755, 7271649, 4.08]
AMPU$$$PERV [280812, 6845115, 4.1]
EPHI$$$VOGE [48604, 7690947, 0.63]
AMPU$$$GRAC [102958, 7271958, 1.42]
GRAC$$$RAJA [121090, 7707372, 1.57]
AMPU$$$RAFF [83766, 7361784, 1.14]
HAMA$$$RAFF [121563, 7782612, 1.56]
EPHI$$$HAMA [102096, 7665465, 1.33]
HEMS$$$PERV [295063, 7307232, 4.04]
DION$$$GRAC [727485, 3528135, 20.62]
JACQ$$$PERV [296332, 7073301, 4.19]
KHAS$$$RAJA [182093, 7632630, 2.39]
GRAC$$$RAFF [106378, 7741056, 1.37]
DION$$$JACQ [714092, 3458850, 20.65]
AMPU$$$HEMS [82186, 7403934, 1.11]
GRAC$$$HEMS [105251, 7779771, 1.35]
RAJA$$$SPGI [116040, 7858041, 1.48]
KHAS$$$RAFF [178798, 7662195, 2.33]
EPHI$$$RAFF [33060, 1047794, 3.16]
PERV$$$RAFF [87329, 978510, 8.92]
GRAC$$$JACQ [34350, 1003792, 3.42]
HAMA$$$VOGE [27679, 1036269, 2.67]
PERV$$$RAJA [86200, 973074, 8.86]
DION$$$HEMS [204787, 460015, 44.52]
HEMS$$$KHAS [52293, 1037710, 5.04]
GRAC$$$KHAS [53878, 1014806, 5.31]
HEMS$$$VOGE [32996, 1055427, 3.13]
KHAS$$$VOGE [52714, 1020056, 5.17]
AMPU$$$DION [188868, 422173, 44.74]
EPHI$$$RAJA [18927, 1044754, 1.81]
JACQ$$$SPGI [33225, 1023105, 3.25]
AMPU$$$SPGI [25895, 998643, 2.59]
KHAS$$$SPGI [53195, 1032263, 5.15]
DION$$$KHAS [202793, 456542, 44.42]
DION$$$RAJA [204758, 459429, 44.57]
EPHI$$$GRAC [34903, 1026073, 3.4]
EPHI$$$HEMS [32659, 1053771, 3.1]
DION$$$EPHI [200823, 449862, 44.64]
JACQ$$$RAFF [32556, 1024401, 3.18]
AMPU$$$RAJA [31127, 993677, 3.13]
HAMA$$$PERV [86910, 967228, 8.99]
HEMS$$$SPGI [9710, 1075554, 0.9]
HAMA$$$KHAS [52778, 1020204, 5.17]
DION$$$VOGE [203450, 456878, 44.53]
SPGI$$$VOGE [34150, 1049736, 3.25]
RAFF$$$SPGI [10053, 1069886, 0.94]
AMPU$$$JACQ [31532, 962652, 3.28]
JACQ$$$RAJA [22552, 1021751, 2.21]
AMPU$$$EPHI [32057, 982012, 3.26]
RAJA$$$VOGE [19144, 1048770, 1.83]
HAMA$$$HEMS [33265, 1056576, 3.15]
EPHI$$$SPGI [33571, 1046038, 3.21]
HAMA$$$JACQ [27442, 1014228, 2.71]
DION$$$RAFF [204699, 459893, 44.51]
JACQ$$$KHAS [51435, 997637, 5.16]
GRAC$$$PERV [87635, 961056, 9.12]
AMPU$$$VOGE [32357, 981443, 3.3]
EPHI$$$KHAS [52190, 1013089, 5.15]
RAFF$$$VOGE [33462, 1050252, 3.19]
GRAC$$$HAMA [35449, 1030218, 3.44]
AMPU$$$KHAS [50329, 966392, 5.21]
AMPU$$$HAMA [32589, 984361, 3.31]
HAMA$$$RAJA [26206, 1046293, 2.5]
AMPU$$$HEMS [24947, 1002698, 2.49]
EPHI$$$PERV [86167, 960065, 8.98]
DION$$$HAMA [202807, 454944, 44.58]
JACQ$$$VOGE [23719, 1014150, 2.34]
PERV$$$VOGE [86515, 964786, 8.97]
EPHI$$$JACQ [23315, 1009793, 2.31]
GRAC$$$VOGE [35360, 1032082, 3.43]
HEMS$$$RAJA [31667, 1067886, 2.97]
RAFF$$$RAJA [32130, 1062258, 3.02]
GRAC$$$SPGI [32633, 1041967, 3.13]
DION$$$PERV [198105, 444591, 44.56]
HAMA$$$SPGI [34192, 1047904, 3.26]
KHAS$$$PERV [82544, 962365, 8.58]
HEMS$$$RAFF [7980, 1081705, 0.74]
DION$$$SPGI [205155, 460768, 44.52]
PERV$$$SPGI [87744, 978054, 8.97]
AMPU$$$PERV [83359, 923113, 9.03]
EPHI$$$VOGE [13017, 1037967, 1.25]
HEMS$$$PERV [87297, 982911, 8.88]
GRAC$$$RAJA [34039, 1039691, 3.27]
AMPU$$$RAFF [25435, 996993, 2.55]
HAMA$$$RAFF [33696, 1050578, 3.21]
EPHI$$$HAMA [27103, 1034031, 2.62]
DION$$$GRAC [201828, 452837, 44.57]
JACQ$$$PERV [85315, 951069, 8.97]
KHAS$$$RAJA [51597, 1028438, 5.02]
GRAC$$$RAFF [31956, 1044871, 3.06]
DION$$$JACQ [197200, 443695, 44.44]
HEMS$$$JACQ [32121, 1028719, 3.12]
GRAC$$$HEMS [31653, 1050241, 3.01]
AMPU$$$GRAC [31079, 983676, 3.16]
RAJA$$$SPGI [32689, 1060681, 3.08]
KHAS$$$RAFF [52519, 1032780, 5.09]

"""

fourfold_degen_codons = ["GCN", "CGN", "GGN", "CTN", "CCN", "TCN", "ACN", "GTN"]
# source: https://en.wikipedia.org/wiki/Codon_degeneracy
# there are 8 4-fold degenerate codons: codons that may contain any base at 1 positions (pos3) and still encode the same amino acid

fourfold_degen_codon_starts = set( [ x.rstrip("N") for x in fourfold_degen_codons ] )

def check_if_fourfold_degenerate (codon_list):
	
	fourfold_degs = [ x for x in codon_list if x[:2] in fourfold_degen_codon_starts ]
#	print fourfold_degs
	if len( fourfold_degs ) == len( [ x for x in codon_list if x != "---"] ): # allows missing codons
		fourfold_deg_column = [ x[2] for x in codon_list ]
#		print "4-fold degenerate"
		return fourfold_deg_column
	else:
#		print "not 4-fold degenerate"
		return None

def read_phylip(INFILE):
	
	indict = {}
	with open(INFILE, "r") as infile:
		infile.readline() # remove header
		for line in infile:
			fields = line.strip("\n").split()
			indict[fields[0]] = fields[1]

	return indict


def write_phylip (concat_dict, outfile):
	
	with open(outfile, "w") as OUTFILE:
		ntaxa = len(concat_dict.keys())
		len_align = len(concat_dict[concat_dict.keys()[0]]) # gets value of first element in dictionary -> this OK since all seqs have same length
		header = str(ntaxa) + " " + str(len_align) + "\n"
		OUTFILE.write(header)
		for sample, seq in concat_dict.items():
			out_line = sample + "    " + seq + "\n"
			OUTFILE.write(out_line)
	OUTFILE.close()


def walk_alignment (indict):
	
	n_columns = len( indict.values()[0] )
	
	tax_order = indict.keys()
	
	fourfold_columns = []
	idx=0
	while idx < len( indict.values()[0] ):
		codon_column = [ indict[x][idx:(idx+3)] for x in tax_order ]
		fourfold_col = check_if_fourfold_degenerate (codon_column)
		fourfold_columns.append( fourfold_col )
#		print fourfold_col
		idx += 3
	
	good_columns = [ x for x in fourfold_columns if x ]
	print len( good_columns )
	
	out_dict = {tax:[] for tax in tax_order}	
	for idx, tax in enumerate( tax_order ):
		for col in good_columns:
			out_dict[tax].append( col[idx] )
	
	
	out_dict_clean = {k: "".join(v) for k,v in out_dict.items() }	
	return out_dict_clean

def count_pairwise_naive_distance ( align_dict ):
	
	all_pairs = set()
	for p in align_dict.keys():
		for q in align_dict.keys():
			if p != q:
				all_pairs.add( "$$$".join( sorted([p,q]) ) )
	
	pairs_dist = {}
	for pair in all_pairs:
		diff = 0
		total = 0
		tax1 = pair.split("$$$")[0]
		tax2 = pair.split("$$$")[1]
		for idx in range(len( align_dict.values()[0] )):
			p = align_dict[tax1][idx]
			q = align_dict[tax2][idx]
			if p != "-" and q != "-":
				total += 1
				if p != q:
					diff += 1
		pairs_dist[pair] = [diff, total, round(float(diff)/float(total)*100, 2) ]
		print pair, pairs_dist[pair]

							
	
	
###############################
# main

infilename = sys.argv[1]

in_alignment = read_phylip(infilename)
#print in_alignment
out_dict_clean = walk_alignment( in_alignment )
count_pairwise_naive_distance ( in_alignment )
count_pairwise_naive_distance ( out_dict_clean )
write_phylip (out_dict_clean, infilename + ".fourfold_deg_sites.txt")

print "Done!"


