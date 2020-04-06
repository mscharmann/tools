from ete2 import Tree
import ete2
import random


def make_nonredundant_pairs (samples):
	
	iterlist1 = samples
	iterlist2 = iterlist1
	
	combinations = set([])
	for a in iterlist1:
		for b in iterlist2:
			if a != b:
				comb = "$".join( sorted([a, b]) )
				combinations.add(comb)
	
#	unique_combs = [list(x) for x in set(tuple(x) for x in combinations)]

	return [ x.split("$") for x in combinations ]

def make_nonredundant_quartets (samples):
	
	## this is a very stupid idea for larger lists of samples!!
	
	iterlist1 = samples
	iterlist2 = samples
	iterlist3 = samples
	iterlist4 = samples
	
	combinations = []
	for a in iterlist1:
		for b in iterlist2:
			if a != b:
				for c in iterlist3:
					if c != a:
						if c != b:
							for d in iterlist4:
								if d != a:
									if d != b:
										if d != c:
											comb = sorted([a, b]) + sorted([c, d])
											combinations.append(comb)
	
	unique_combs = [list(x) for x in set(tuple(x) for x in combinations)]

	return unique_combs

def extract_contrast_pairs_from_quartet (q):
	
	X, Y, Z, W = q[0], q[1], q[2], q[3]
	
	return set( [ "$".join(sorted( [X, Z])), "$".join(sorted( [X, W] )), "$".join(sorted( [Y, Z] )), "$".join(sorted( [Y, W] )) ] )
	

def get_tip_names ( t_node ):
	
	leaves = t_node.get_leaves(is_leaf_fn=None)
	tips = []
	for l in leaves:
		tips.append( l.name)
	
	return tips

def cull_down_quartets ( good_quartets, max_n_contrasts ):
	
	"""
	each quartet X,Y,Z,W contains the four pair-comparisons X-Z, X-W, Y-Z, Y-W
	
	"""
	
	contrasts_dict = {}
	for q in good_quartets:
		X, Y, Z, W = q[0], q[1], q[2], q[3]
		try:
			contrasts_dict[ "$".join(sorted( [X, Z] )) ].append( q )
		except KeyError:
			contrasts_dict[ "$".join(sorted( [X, Z] )) ] = [ q ]
		try:
			contrasts_dict[ "$".join(sorted( [X, W] )) ].append( q )
		except KeyError:
			contrasts_dict[ "$".join(sorted( [X, W] )) ] = [ q ]
		try:
			contrasts_dict[ "$".join(sorted( [Y, Z] )) ].append( q )
		except KeyError:
			contrasts_dict[ "$".join(sorted( [Y, Z] )) ] = [ q ]
		try:
			contrasts_dict[ "$".join(sorted( [Y, W] )) ].append( q )
		except KeyError:
			contrasts_dict[ "$".join(sorted( [Y, W] )) ] = [ q ]
	
	# the downsampling:
	culled_quartets = set([])
	for pair in contrasts_dict.keys():
		try:
			accepted_contrasts = random.sample(contrasts_dict[pair], max_n_contrasts)
		except ValueError:
			accepted_contrasts = contrasts_dict[pair]
		for a in accepted_contrasts:
			culled_quartets.add( "$".join(a) )
	
	print culled_quartets
	
	final_culled_quartets = [x.split("$") for x in culled_quartets]
	
	print """after downsampling to max. {0} contrasts per species pair,
retained {1} ABBA-BABA quartets ( {2} % )""".format( max_n_contrasts , len( final_culled_quartets ) , float(len(final_culled_quartets))/float(len( good_quartets ))*100 )
		
	return final_culled_quartets


def sort_and_output ( quartets ):
	
	## sort the outcome and write to file::
	quart_dict = {}
	for i in  quartets:
		quart_dict["".join(i)] = i
	
	outlines = ["\t".join(["X", "Y", "Z", "W"]) ]
	for i in sorted(quart_dict.keys()):
		outlines.append( "\t".join( quart_dict[ i ] ))
	
	with open("ABBA-BABA_quartets.ExaML_result.T1.1.txt", "w") as OUTFILE:
		OUTFILE.write("\n".join(outlines))


def check_if_conforming_to_species_tree ( quartets , species_tree):
	
	# Is the condition ((p1,p2),(p3,p4)) true on the species tree?
	
	good_quartets = []
	for q in quartets:
		p1, p2, p3, p4 = q[0], q[1], q[2], q[3]
		if len( set([p1,p2,p3,p4] )) == 4:
			tree_hypothesis = Tree( "((" + p1 + "," + p2 + "),(" + p3 + "," + p4 + "));" )
			tree_comps = species_tree.robinson_foulds(tree_hypothesis,unrooted_trees=True)
			robinson_foulds_dist = tree_comps[0]
			if robinson_foulds_dist == 0: ## means that tree_hypothesis is a subtree of the species tree!
				good_quartets.append([p1,p2,p3,p4])
		else:
			print [p1,p2,p3,p4]
				
	return good_quartets


def make_quartets_treewalking (t):

	"""	
		# idea: starting from outgr, get next node, get all 
		for internal_node in allnodes:
			if len(leaves_from_node) >= 3:
				left_leaves =
				right_leaves =
				all_left_pairs = nonredundant_pairs
				all_right_pairs = nonredundant_pairs
				left_triplets = all_left_pairs * right_leaves
				right_triplets = all_right_pairs * left_leaves
			
				quartets.append( [  ] ) 
	"""
	all_tips = get_tip_names( t )
	
	
	
	## we search more constrasts until all possible_nonredundant_pairs have been found among them:
	all_valid_quartets = set( )
	all_contrast_pairs = set()
	
	# stop searching for more pairs when the number didn't change for 50 cycles:
	nonr_pairs_count_history = range(100) ## fill recorder with nonsense
	
	while len(set( nonr_pairs_count_history[-100:] )) > 1:

		nonr_pairs_count_history.append( len( all_contrast_pairs ) )

		a = random.sample(all_tips, 1) ## choose a random root
		outgr = a[0]
		t.set_outgroup( outgr )
		print len( all_contrast_pairs )
		b = random.sample(t.get_descendants(), len(all_tips)/2 ) # chose random internal (and terminal = leaf) nodes
		for node in b: # iterate over  nodes
			if len( node.get_descendants() ) >= 3: # reject leaves and nodes that have too few leaves on them
				# distinguish the two descendant subtrees:
				left_subtree = node.children[0]
 				right_subtree = node.children[1]
			
				left_tips = get_tip_names( left_subtree )
				right_tips = get_tip_names( right_subtree )
			
				# get all possible leaf2, leaf1 pairings:
				all_left_pairs = make_nonredundant_pairs ( left_tips )
				all_right_pairs = make_nonredundant_pairs ( right_tips )
				
				try:
					subsampled_left = random.sample(all_left_pairs, 1+(len(all_tips)/20) )
				except ValueError:
					subsampled_left = all_left_pairs
				try:
					subsampled_right = random.sample(all_right_pairs, 1+(len(all_tips)/20) )
				except ValueError:
					subsampled_right = all_right_pairs
				
				for l in subsampled_left:
					for r in subsampled_right:
#						print l, r, "$".join( sorted(["$".join(l), "$".join(r)]) ) 
#						quartetlist.append( "$".join( l + r  ))
						cps = extract_contrast_pairs_from_quartet ( l + r )
						if not cps.issubset(all_contrast_pairs):
							all_contrast_pairs = all_contrast_pairs.union( cps )
							all_valid_quartets.add( "$".join( sorted(["$".join(l), "$".join(r)])  ) )

				
#				all_valid_quartets = all_valid_quartets.union(set(quartetlist))
				
#				# get the triplets:
#				quartetlist = []
#				for l in all_left_pairs:
#					for r in right_tips:
#						all_valid_quartets.add( "$".join( sorted([outgr, r]) + l )  )
#						quartetlist.append( "$".join( sorted([outgr, r]) + l  ))
#				all_valid_quartets = all_valid_quartets.union(set(quartetlist))
			
#				right_triplets = []
#				for r in all_right_pairs:
#					for l in left_tips:
#						all_valid_quartets.add( "$".join( sorted([outgr, l]) + r ) )
			
# 			print left_subtree
# 			print left_triplets
# 			print right_subtree
# 			print right_triplets
# 			print "-----------------------------"
	
#	print all_valid_quartets
	
	print "found {0} valid ABBA-BABA quartets".format( len( all_valid_quartets ) )
	
	return all_valid_quartets


			
############# MAIN

treefile = "biolmiss_12sites.25.0.1.rexptrunc_0.33_genos.completeness_0.6.ML_tree.txt"
#treefile = "SNPs.present_0.5perpop.minpop2.poplevel.majority_consensus.SH_support.figtree.txt"
#treefile = "toytree.txt"
treefile="ExaML_result.T1"
max_n_contrasts = 1

t = Tree( treefile )
t.unroot()

print t

quartets_from_heuristic = make_quartets_treewalking (t)

quartets_from_heuristic = [x.split("$") for x in quartets_from_heuristic]

sort_and_output ( quartets_from_heuristic )

exit()

all_tip_names = get_tip_names( t )
#all_pairs_of_tips = make_nonredundant_pairs ( all_tip_names )
#joined_pairs = ["@".join(x) for x in all_pairs_of_tips]
#print joined_pairs
#print len( joined_pairs )
#raw_quartets = make_nonredundant_pairs ( joined_pairs )
#
#print "raw quartets:	", len( raw_quartets )
#
#cleaned_quartets = []
#for q in raw_quartets:
#	p1, p2 = q[0].split("@")	
#	p3, p4 = q[1].split("@")
#	if len( set([p1,p2,p3,p4] )) == 4:
#		cleaned_quartets.append( [p1,p2,p3,p4] )

cleaned_quartets = make_nonredundant_quartets (all_tip_names)	
	
print "cleaned quartets:	", len( cleaned_quartets )

good_quartets = check_if_conforming_to_species_tree ( cleaned_quartets, t )

print "good quartets:	", len( good_quartets )

## Now, we want that each "pairwise" comparison is conducted only max_n_contrasts times.
culled_quartets = cull_down_quartets ( good_quartets, max_n_contrasts )


sort_and_output ( culled_quartets )

print "Done!"
