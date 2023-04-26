import pandas as pd
from itertools import combinations
from collections import Counter
import numpy as np

tissues_fname = 'saved_gtex_datasets.txt'
corr_thres = 0.6

# load tissues types for which module files will be available
tissues = [ l.strip().split('\t')[1].strip('"') for l in open(tissues_fname).readlines()[1:] ]

data = {}

# 1-19 since we have 19 tissues
min_overlaps = [10]

for tissue in tissues:
	fname = "corrs_%s.txt" % tissue
	print "Loading", tissue, " ... ",

	# load correlation matrix. We load one at the time to avoid
	# potential memory limitations. These matrices are big - 
	# they take approx 1.6GB each on my system
	corrs = pd.read_csv(fname, sep='\t', header=None, skiprows=1).astype(float).values
	np.fill_diagonal(corrs, 0.0)

	# load labels
	labels = open(fname).readline().strip().split("\t")

	# check there are no duplicated labels
	noempty_labels = open(fname).readline().strip().split()
	assert(len(noempty_labels) == len(list(set(noempty_labels))))
	
	# get matrices indices with absolute correlation values higher than threshold
	pair_idxs = zip(*np.where(np.abs(corrs) >= corr_thres))

	# get corresponding label pairs - discarding cases in which the label is empty
	pairs = [ frozenset((labels[i[0]], labels[i[1]])) for i in pair_idxs if labels[i[0]] != '' and labels[i[1]] != '']
	
	# eliminate duplicates
	data[tissue] = list(set(pairs))

# all_interactions will contain all the pairs for all the tissues - with duplicates if they are found
# in more than one tissue. Allows us to calculate in how many tissues a given pair is present, even
# though we lose the pair-tissue association

all_interactions = []
for v in data.itervalues():
	all_interactions.extend(v)

print "all interactions with duplicates:", len(all_interactions)
print "all unique interactions:", len(set(all_interactions))

# for each interaction, count in how many tissue it is present
counted_interactions = Counter(all_interactions)

# count how many times we get 1, 2, 3, 4 ... n tissues per interaction
counted_histogram = Counter(counted_interactions.values())

# prepare the cumulative sum
cumulative = []

for i in range(1, len(tissues)+1):
	cumulative.append(counted_histogram[i])

for i,n in enumerate(cumulative):
	print "%d\t%d" % (i+1, sum(cumulative[i:]))

# write one file per minimum number of tissue. Each file will be a IID-like
# file with pairs of interactions, each of them with a different minimum number
# of tissues in which said interaction is found

for min_overlap in min_overlaps:
	filtered_dataset = []

	# generate specific dataset filtering for minimum number of tissues in which 
	# interactions are present 
	for i,v in counted_interactions.iteritems():
		if v >= min_overlap:
			filtered_dataset.append(i)
	filtered_dataset = zip(*filtered_dataset)

	if len(filtered_dataset) == 0:
		print "no data with minimum %d interaction; skipping" % min_overlap
		continue

	# generate dummy data to keep the same structure as IID
	up1 = ['xxxxxx' for i in filtered_dataset[0]]
	up2 = up1
	methods = ['coexpression' for i in filtered_dataset[0]]
	pmids = ['' for i in filtered_dataset[0]]
	dbs = pmids
	evidence_type = ['exp' for i in filtered_dataset[0]]

	# create list with number of tissues for each interaction
	ntissues = [ counted_interactions[frozenset([filtered_dataset[0][i], filtered_dataset[1][i]])] for i in range(len(filtered_dataset[0])) ]

	# create dataframe and write it to disk as csv file
	df_filtered_dataset = pd.DataFrame(	{'uniprot1' : up1,
						'uniprot2'  : up2,
						'symbol1'   : filtered_dataset[0],
						'symbol2'   : filtered_dataset[1],
						'methods'   : methods,
						'pmids'     : pmids,
						'dbs'       : dbs,
						'evidence type' : evidence_type,
						'ntissues' : ntissues})

	df_filtered_dataset.to_csv('coexpression-interactions_%.1f_min%s.csv' % (corr_thres, min_overlap))	
	
