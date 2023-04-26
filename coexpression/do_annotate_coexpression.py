# do-db.py - filter db database according to selected criteria
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
import sys


db_fname = 'coexpression-interactions_0.6_min10.csv'
go_association_fname =  'goa_human.gaf'
go_fname = 'go.obo'

go_dna_repair = 'GO:0006281'
go_response_dna_damage = 'GO:0006974'

def parse_db(fname):
	"""parse the database file and
	1. remove entries in which the two genes are the same
	"""

	d = pd.read_table(fname, sep=',')
	d.update(d['methods'].str.split(';'))
	d.update(d['evidence type'].str.split(';'))
	d = d[ d.apply(lambda x: x['symbol1'] != x['symbol2'], axis=1)  ]
	d = d.set_index(['symbol1', 'symbol2'])
	return d

def parse_go(fname, children_go_dna_repair, children_go_response_dna_damage):
	"""parse GO .goa file with association between GO and proteins. 
           Add the DNA repair or response DNA damage go if children of those two
           are present.
        """

        objid = 'DB_Object_Symbol'
        goid  = 'GO_ID'
        syno  = 'Synonym'
        data  = defaultdict(set)

        with open(fname) as fh:
            proteins = GOA.gafiterator(fh)
            for p in proteins:
                data[p[objid]].add(p[goid])
                for s in p[syno]:
                	data[s].add(p[goid])
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

# Parse db filtering out prediction PMIDs
db = parse_db(db_fname)

# Parse GO annotation for proteins
gos = parse_go(go_association_fname, children_go_dna_repair, children_go_response_dna_damage)

# Add GO annotation to db
db = add_go(db, gos)

# Filter db by GO terms
db = db[ (db.apply(lambda x:      (go_dna_repair in x['go1'] and not (go_dna_repair in x['go2'] or go_response_dna_damage in x['go2']) ) \
                                  or (go_dna_repair in x['go2'] and not (go_dna_repair in x['go1'] or go_response_dna_damage in x['go1']) ), 
                                  axis=1 ))]

# Change the format of the database to make it more human-friendly
print(db)

db.update(db['methods'].str.join(', '))
#db.update(db['pmids'].str.join(', '))
#db.update(db['dbs'].str.join(', '))
db.update(db['evidence type'].str.join(', '))
db.update(db['go1'].str.join(', '))
db.update(db['go2'].str.join(', '))

# Save the final database as csv
db.to_csv(db_fname + '.go_annotated.csv')
