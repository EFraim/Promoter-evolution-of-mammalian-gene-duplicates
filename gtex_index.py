#!env python3
from pandas import read_csv
import shelve
s=shelve.open('/data/db/import/gtex_dict_pilot.db', 'n')
for chunk in read_csv('/data/db/import/GTEx_Analysis_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm__Pilot_V3_patch1.gct', sep='\t', header=2, chunksize=10000):
  for l in chunk.itertuples():
    s[l[1].split('.')[0]]=l[3:]

s.sync()
s.close()

