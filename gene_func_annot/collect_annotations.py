
import gzip, sys

"""
DONE GO terms TAIR
DONE GO terms SwissProt

DONE InterPro domains
DONE GO terms Interpro

DONE human-readable TAIR
DONE human-readable SwissProt
DONE human-readable Interpro

# before this works, download:

wget http://current.geneontology.org/annotations/tair.gaf.gz

wget http://current.geneontology.org/annotations/goa_uniprot_all.gaf.gz # this is more complete, but takes very long to search; is 16 GB zipped..
wget http://current.geneontology.org/annotations/goa_uniprot_all_noiea.gaf.gz # this is less complete, but VERY quick, only 6 Mbp zipped.


wget https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/interpro2go

# PlantTFDB v5.0 Plant Transcription Factor Database
wget http://planttfdb.gao-lab.org/download/TF_list/Ath_TF_list.txt.gz


"""




def parse_all_protein_IDs (infasta):
	
	print("parsing input fasta ", infasta)
	
	all_IDs = set()
	with open(infasta, "r") as I:
		for line in I:
			if line.startswith(">"):
				id = line.lstrip(">").strip()
				all_IDs.add(id)
	return(sorted(list(all_IDs)))


def parse_TF_database (infile):
	
	# besthit by evalue.
	
	print("parsing transcription factor database ", infile)
	
	ATids_to_TFs = {}
	with gzip.open(infile, "rt") as I:
		for line in I:
			if len(line) > 1:
				fields = line.strip("\n").split("\t")
				ATids_to_TFs[fields[1]] = fields[2]
	
	return(ATids_to_TFs)


def parse_blast_besthit (infile):
	
	# besthit by evalue.
	
	print("parsing BLAST hits ", infile)
	
	blast_besthits = {}
	with gzip.open(infile, "rt") as I:
		for line in I:
			if len(line) > 1:
				fields = line.strip().split("\t")
				qid = fields[0]
				hit_id = fields[1]
				evalue = float(fields[10])
				bitscore = float(fields[11])
				if evalue <= 1e-5:
					try:
						prev_bitscore = blast_besthits[qid][1]
						if prev_bitscore < bitscore:
							blast_besthits[qid] = [hit_id, bitscore]
					except KeyError:
						blast_besthits[qid] = [hit_id, bitscore]
	
	return(blast_besthits)	


def	parse_TAIR_GAF_for_human_readable_description_and_GO_terms (gaf_infile):
	
	# http://geneontology.org/docs/go-annotation-file-gaf-format-2.1/
	
	print("parsing TAIR GAF ", gaf_infile)
	
	ids_and_description = {}
	ids_and_GOs = {}
	with gzip.open(gaf_infile, "rt") as I:
		for line in I:
			if not line.startswith("!"):
				fields = line.strip().split("\t")
				ids_and_description[ fields[1] ] = "|".join( fields[10].split("|")[1:] )
				try:
					ids_and_GOs[ fields[1] ].append(fields[4])
				except KeyError:
					ids_and_GOs[ fields[1] ] = [fields[4]]
	
	tmp = {}
	for k,v in ids_and_GOs.items():
		tmp[k] = ",".join(list(set(v)))
	
	
	return(ids_and_description,tmp)
		

def	parse_uniprot_GAF_for_human_readable_description_and_GO_terms (gaf_infile):
	
	# http://geneontology.org/docs/go-annotation-file-gaf-format-2.1/
	
	print("parsing Uniprot GAF ", gaf_infile)
	
	ids_and_description = {}
	ids_and_GOs = {}
	with gzip.open(gaf_infile, "rt") as I:
		for line in I:
			if not line.startswith("!"):
				fields = line.strip().split("\t")
				ids_and_description[ fields[1] ] = fields[9]
				try:
					ids_and_GOs[ fields[1] ].append(fields[4])
				except KeyError:
					ids_and_GOs[ fields[1] ] = [fields[4]]
	
	tmp = {}
	for k,v in ids_and_GOs.items():
		tmp[k] = ",".join(list(set(v)))
	
	
	return(ids_and_description,tmp)


def parse_interproscan_hits (infile):
	
	print("parsing InterProscan hits ", infile)
	
	ipr_hits = {}
	with gzip.open(infile, "rt") as I:
		for line in I:
			if len(line) > 1:
				fields = line.strip().split("\t")
				qid = fields[0]
				ipr_hit = fields[11]
				descr = fields[12]
				if ipr_hit != "-":
					try:
						ipr_hits[qid].append( ipr_hit + " (" + descr + ")" ) 
					except KeyError:
						ipr_hits[qid] = [ ipr_hit + " (" + descr + ")" ]

	tmp = {}
	for k,v in ipr_hits.items():
		tmp[k] = "; ".join(list(set(v)))

	
	return(tmp)	


def parse_interpro_GOs (infile):
	
	print("parsing InterPro GOs ", infile)
	
	ids_and_GOs = {}
	with open(infile, "r") as I:
		for line in I:
			if not line.startswith("!"):
				id = line.split()[0].split(":")[1]
				go = line.strip().split(" ; ")[-1]
				# print(id,go)
				try:
					ids_and_GOs[ id ].append(go)
				except KeyError:
					ids_and_GOs[ id ] = [go]
	
	tmp = {}
	for k,v in ids_and_GOs.items():
		tmp[k] = ",".join(list(set(v)))
	
	
	return(tmp)


def parse_GTF (infile):
	
	print("parsing GTF ", infile)
	
	genes_to_coords_dict = {}
	with open(infile, "r") as I:
		for line in I:
			if not line.startswith("#"):
				fields = line.strip("\n").split("\t")
				if fields[2] == "transcript":
					c = fields[0]
					s = fields[3]
					e = fields[4]
					gene = fields[8].split()[3].strip('"')
					genes_to_coords_dict[gene] = "\t".join([c,s,e])

	return(genes_to_coords_dict)

### MAIN



input_seqids = parse_all_protein_IDs (sys.argv[1])

genes_to_coords_dict = parse_GTF ("augustus.hints.gtf")

ATids_to_TFs = parse_TF_database ("Ath_TF_list.txt.gz")

sprot_besthits = parse_blast_besthit ("blastresult.SProt.txt.gz")
TAIR_besthits = parse_blast_besthit ("blastresult.TAIR.txt.gz")


tair_descriptions, tair_GOs = parse_TAIR_GAF_for_human_readable_description_and_GO_terms ("tair.gaf.gz")
sprot_descriptions, sprot_GOs = parse_uniprot_GAF_for_human_readable_description_and_GO_terms ("goa_uniprot_all_noiea.gaf.gz")

ipr_hits_with_descr = parse_interproscan_hits ("interproscan_result.txt.gz")
# print(ipr_hits_with_descr)

ipr_GOs = parse_interpro_GOs ("interpro2go")
#print(ipr_GOs)

print("collecting output")

outlines = []
outlines.append("\t".join(["query_ID","chr","start","end","BLAST_hit_TAIR","TAIR_description","AT_transcription_factor_family","BLAST_hit_SwissProt","SwissProt_description","InterProscan_hits_(Description)","TAIR_GOterms","SwissProt_GOterms","InterProscan_GOs"]))
for s in input_seqids:
	coords = genes_to_coords_dict[ ".".join(s.split(".")[:-1]) ]
	try:
		tairhit = TAIR_besthits[s][0].split(".")[0]
	except KeyError:
		tairhit = "NA"
	try:
		sprot_hit = sprot_besthits[s][0].split("|")[1]
	except KeyError:
		sprot_hit = "NA"
	try:
		tair_descr = tair_descriptions[tairhit]
	except KeyError:
		tair_descr = "NA"
	try:
		tgs = tair_GOs[tairhit]
	except KeyError:
		tgs = ""
	try:
		sprot_descr = sprot_descriptions[sprot_hit]
	except KeyError:
		sprot_descr = "NA"
	try:
		spgs = sprot_GOs[sprot_hit]
	except KeyError:
		spgs = ""

	try:
		ipr_hits = ipr_hits_with_descr[s]
	except KeyError:
		ipr_hits = "NA"

	if ipr_hits != "NA":
		xa = ipr_hits.split("; ")
		xb = [x.split()[0] for x in xa]
		ipr_gos = []
		for l in xb:
			#print(l)
			try:
				ipr_gos.append(ipr_GOs[l])
			except KeyError:
				None
		if len(ipr_gos) > 0:
			iprgos_out = ",".join( list(set(",".join(ipr_gos).split(","))) )
		else:
			iprgos_out = ""	 		
	else:
		iprgos_out = ""
	
	if tairhit != "NA":
		try:
			TF_fam = ATids_to_TFs[tairhit]
		except KeyError:
			TF_fam = ""
	else:
		TF_fam = ""
		
	
	outline = "\t".join([s,coords,tairhit,tair_descr,TF_fam,sprot_hit,sprot_descr,ipr_hits,tgs,spgs,iprgos_out])
	outlines.append(outline)
	
	

with open(sys.argv[1] + ".functional_annotations.txt", "w") as O:
	O.write( "\n".join(outlines) + "\n" )

print ("Done!")

	
	
