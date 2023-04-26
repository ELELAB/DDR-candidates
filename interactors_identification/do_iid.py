# do-iid.py - filter iid database according to selected criteria
#
#     Copyright (C) 2018 Matteo Tiberti <matteo.tiberti@gmail.com>
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

from collections import defaultdict
import pandas as pd
from Bio.UniProt import GOA
import time
from goatools import obo_parser


iid_fname = 'iid.human.2018-05.txt'
go_association_fname =  'goa_human.gaf'
go_fname = 'go.obo'

go_dna_repair = 'GO:0006281'
go_response_dna_damage = 'GO:0006974'
prediction_pmids = set(['25402006', '23023127', '21836163', '16082366', 'doi:10.1007/978-3-540-73060-6_4'])

def parse_iid(fname, prediction_pmids):
	"""parse IID database file and remove
		1. entries in which proteins interact with themselves
		2. PMIDs relative to prediction studies (as we won't need them)
	"""

	d = pd.read_table(fname, usecols=range(0,8), quoting=3)
	d.update(d['methods'].str.split(','))
	d.update(d['pmids'].str.split(','))
	d.update(d['dbs'].str.split(','))
	d.update(d['evidence type'].str.split(','))
        pmids_filt = {'pmids':[]}
	for entry in d['pmids']:
                pmids_filt['pmids'].append( sorted(list(set(entry) - prediction_pmids)))
        d.update(pd.DataFrame(pmids_filt))
	d = d[ d.apply(lambda x: x['uniprot1'] != x['uniprot2'], axis=1)  ]
	d = d.set_index(['uniprot1', 'uniprot2'])
	return d

def parse_go(fname, children_go_dna_repair, children_go_response_dna_damage):
	"""parse GO .goa file with association between GO and proteins. 
           Add the DNA repair or response DNA damage go if children of those two
           are present.
        """

        objid = 'DB_Object_ID'
        goid  = 'GO_ID'
        syno  = 'Synonym'
        data  = defaultdict(set)

        with open(fname) as fh:
                proteins = GOA.gafiterator(fh)
                for p in proteins:
                        data[p[objid]].add(p[goid])
	for k,v in data.iteritems():
		if not v.isdisjoint(children_go_dna_repair):
			data[k].add(go_dna_repair)
		if not v.isdisjoint(children_go_response_dna_damage):
			data[k].add(go_response_dna_damage)
        return data

def add_go(data, gos):
	"""add GO annotation to the interaction database"""
	go1 = []
	go2 = []
	for idx, row in data.iterrows():
		go1.append(gos[idx[0]])
		go2.append(gos[idx[1]])
	data = data.assign(go1=pd.Series(go1).values)
	data = data.assign(go2=pd.Series(go2).values)
	return data

# Get children of key GO terms from GO
gos = obo_parser.GODag(go_fname, optional_attrs=['relationships'])

children_go_dna_repair = set([go_dna_repair]) | gos[go_dna_repair].get_all_children()
children_go_response_dna_damage = set([go_response_dna_damage]) | gos[go_response_dna_damage].get_all_children()

# Parse IID filtering out prediction PMIDs
iid = parse_iid(iid_fname, prediction_pmids)

# Parse GO annotation for proteins
gos = parse_go(go_association_fname, children_go_dna_repair, children_go_response_dna_damage)

# Add GO annotation to IID
iid = add_go(iid, gos)

# Filter IID by GO terms
iid = iid[ (iid.apply(lambda x:     (go_dna_repair in x['go1'] and not (go_dna_repair in x['go2'] or go_response_dna_damage in x['go2']) ) \
                                  or (go_dna_repair in x['go2'] and not (go_dna_repair in x['go1'] or go_response_dna_damage in x['go1']) ), 
                                  axis=1 ))]

# Filter IID by keeping only interactions having experimental evidence
iid = iid[ iid.apply(lambda x: 'exp' in x['evidence type'], axis=1) ]

# Add column with the number of PMIDs to the dataframe - notice that we 
# already remove the PMIDs relative to the prediction studies when
# parsing the database
iid.insert(4, "n_pmids", iid.apply(lambda x: len(x['pmids']), axis=1))

# Filter out entries with just 1 PMID
iid = iid[ iid.apply(lambda x: x['n_pmids']  >= 2, axis=1) ]

# Sort the entries in the database by number of PMIDs
iid = iid.sort_values(by='n_pmids', ascending=False)

# Change the format of the database to make it more human-friendly
iid.update(iid['methods'].str.join(', '))
iid.update(iid['pmids'].str.join(', '))
iid.update(iid['dbs'].str.join(', '))
iid.update(iid['evidence type'].str.join(', '))
iid.update(iid['go1'].str.join(', '))
iid.update(iid['go2'].str.join(', '))

# Save final database as csv file
iid.to_csv('iid_interactions.csv')

