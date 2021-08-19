from ete3 import Tree
from ete3 import EvolTree
import ete3
import random
import sys
import os
import numpy as np

# implements ete3 EvolTree b_free (the free-branch model: where omega-foreground and omega-background are free and different)
# corresponds to the lysozyme example in pamlDOC
# checks if there are at least 1 synonymous substitutions expected on this branch... otherwise not very meaningful to estimate omega (it will be very high).

def get_tip_names ( t_node ):
	leaves = t_node.get_leaves(is_leaf_fn=None)
	tips = []
	for l in leaves:
		tips.append( l.name)	
	return tips
	

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


def make_nonredundant_pairs (samples):
	iterlist1 = samples
	iterlist2 = iterlist1	
	combinations = []
	for a in iterlist1:
		for b in iterlist2:
			if a != b:
				comb = sorted([a, b])
				combinations.append(comb)
	
	unique_combs = [list(x) for x in set(tuple(x) for x in combinations)]
	return unique_combs


def make_clean_fasta (seqids, seqdatadict):
	outseq_dict = {k:"" for k in seqids}
	inseqs = [seqdatadict[k] for k in seqids]
	snps = 0
	for pos in range(len(seqdatadict[seqids[0]])):
		columnset = set([x[pos] for x in inseqs])	
		if not "-" in columnset:
			for s in seqids:
				outseq_dict[s] += seqdatadict[s][pos]
				if len(columnset) > 1:
					snps += 1
		
	if len(outseq_dict[seqids[0]]) > 30 and snps >= 3:	
		outlines = []
		for k,v in outseq_dict.items():
			outlines.append(">" + k)
			outlines.append(v)
		return "\n".join(outlines)
	else:
		return None


def export_subtree (tree, wanted_tips):
	outfilename = "__".join(wanted_tips) + ".tre"
	tree.prune(wanted_tips, preserve_branch_length = True)
	print(maintree)

############## MAIN

treefile = sys.argv[1]

#treefile = "OG0000209.CDS.aln.tre"

OG_name = treefile.rstrip(".CDS.aln.tre")
fastafile = OG_name + ".CDS.aln"

t = Tree( treefile )
t.unroot()

all_leaf_names = get_tip_names (t)

all_species = sorted(list(set([ x.split("_")[0] for x in all_leaf_names])))

seqdatadict = read_fasta (fastafile)

# extract all pairwise distances between leaves:
all_leaf_pairs = make_nonredundant_pairs(all_leaf_names)

sp_omegas = []
for sp in all_species:
	print(sp)
	omega_list = []
	list_of_tempdirs = []
	# get all seqs of that species
	seqids_of_species = [x for x in all_leaf_names if x.split("_")[0] == sp ]
	seqids_of_other_species = [x for x in all_leaf_names if x.split("_")[0] != sp ]
	if len(seqids_of_other_species) < 3:
		sp_omegas.append("NA")
		print ("not enough other species sequences, skipping ", sp)
		continue
	# find the three closest sequences in the gene tree that belong to different species: use topology only
	for seqid in seqids_of_species:
		all_dists = []
		for othersp_seq in seqids_of_other_species:
			dist = (t&seqid).get_distance(othersp_seq, topology_only=True) 
			all_dists.append( dist )
		# find indexes of the three shortest distances
		try:
			idxes_of_3_smallest = np.argpartition( np.array( all_dists ) , 3)[:3]
		except ValueError:
			idxes_of_3_smallest = np.argpartition( np.array( all_dists ) , 2) # for the case that list is only 3 items long
		closest_seq_ids = [seqid]
		for d in idxes_of_3_smallest:
			closest_seq_ids.append( seqids_of_other_species[d])
		# ete3 has codeml handling implemented!! No need for own functions.
		subtree = t.copy()
		subtree.prune(closest_seq_ids, preserve_branch_length = True)
		subtree.unroot()
		evotree = EvolTree( subtree.write() )
		subfasta = make_clean_fasta(closest_seq_ids, seqdatadict)
		if not subfasta:
			omega_list.append( "NA" )
			continue
		else:
			evotree.link_to_alignment (subfasta)			
			workdirname = './codeml_' + "__".join(closest_seq_ids)
			evotree.workdir = workdirname
			list_of_tempdirs.append( workdirname ) 
			# mark the foreground branch
			foreground_leafnode = evotree&seqid
#			print (seqid)
#			print(foreground_leafnode.node_id)
#			print (evotree.write())
			evotree.mark_tree ([foreground_leafnode.node_id], ['#1'])
#			print (evotree.write())
							
			evotree.run_model ('b_free.run')
			b_free_fit = evotree.get_evol_model('b_free.run')
			out_branches_dict = b_free_fit.branches
			for b in out_branches_dict:
				if out_branches_dict[b]["mark"] == " #1":
					# check if there are at least 1 synonymous substitutions expected on this branch... otherwise not very meaningful to estimate omega (it will be very high).
					if out_branches_dict[b]["S"] * out_branches_dict[b]["dS"] >= 1.0: 
						omega = out_branches_dict[b]["w"]
					else:
						omega = "NA"
					break
			omega_list.append( omega )
	numeric_omegas = [float(x) for x in omega_list if not x == "NA"]
	try:
		avg_omega = sum(numeric_omegas)/float(len(numeric_omegas))
	except ZeroDivisionError:
		avg_omega = "NA"
	sp_omegas.append(avg_omega)
	# cleanup
	cmd = "rm -r"
	for dir in list_of_tempdirs:
		cmd += " " + dir 
	os.system(cmd)

with open(treefile + ".rep.txt", "w") as OUTF:
	outlines = []
	for p,q in zip(all_species,sp_omegas):
		outlines.append(p + "\t" + str(q))
	OUTF.write("\n".join(outlines)+"\n")
		