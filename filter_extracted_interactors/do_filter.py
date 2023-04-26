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

def filter_bygene(x, genelist):
    gene_name = None
    if go_dna_repair in x['go1'].split(", "):
    	gene_name = x['symbol2']
    if go_dna_repair in x['go2'].split(", "):
        gene_name = x['symbol1']
    if gene_name == None:
        print(x)
        exit(1)
    if gene_name in genelist:
        return False
    return True

from collections import defaultdict
import pandas as pd
from Bio.UniProt import GOA
import time
from goatools import obo_parser

go_dna_repair = 'GO:0006281'
go_response_dna_damage = 'GO:0006974'

interaction_files = ['extracted_interactors.csv']
crapome_file = 'crapome_matrix_reduced_freq.csv'
crapome = pd.read_csv('crapome_matrix_reduced_freq.csv', sep=';')

housekeeping = [ i.strip() for i in open('housekeeping_genes.txt') ]

crapome_all  = crapome['GENEID'].values
crapome_10   = crapome[ crapome['FREQ'] >= 0.1 ]['GENEID'].values
crapome_50   = crapome[ crapome['FREQ'] >= 0.5 ]['GENEID'].values

print(len(housekeeping))
print(len(crapome_all))
print(len(crapome_10))
print(len(crapome_50))

print('fname\tfull\thk\thk+C50\thk+C10\thk+Call')
for fname in interaction_files:

	try:
		iid = pd.read_csv(fname)
	except IOError:
		continue

	print(fname),
	
	iid = iid.fillna('')
	print(len(iid), '\t'),

	iid = iid[ iid.apply(filter_bygene, axis=1, genelist=housekeeping) ]
	iid.to_csv(fname[:-4] + '_housekeeping' + '.csv')
	print(len(iid), '\t'),
	
	#print iid

	iid2 = iid[ iid.apply(filter_bygene, axis=1, genelist=crapome_50)  ]
	iid2.to_csv(fname[:-4] + '_housekeeping_crapome_50' + '.csv')
	print(len(iid2)), '\t',

	#print iid2

	iid2 = iid[ iid.apply(filter_bygene, axis=1, genelist=crapome_10)  ]
	iid2.to_csv(fname[:-4] + '_housekeeping_crapome_10' + '.csv')
	print(len(iid2)), '\t',

	#print iid2

	iid2 = iid[ iid.apply(filter_bygene, axis=1, genelist=crapome_all) ]
	iid2.to_csv(fname[:-4] + '_housekeeping_crapome_all' + '.csv')
	print(len(iid2)), '\t'

	#print iid2

