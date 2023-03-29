from gffutils import FeatureDB
from sys import argv
from itertools import chain

db=FeatureDB(argv[1])
for t in chain(db.features_of_type('mRNA'), db.features_of_type('pseudogenic_transcript'), db.features_of_type('V_gene_segment'), db.features_of_type('J_gene_segment'), db.features_of_type('C_gene_segment')):
    print("%s\t%s"%( t.id.split(':')[-1]+'.'+t.attributes['version'][0],t.attributes['Parent'][0].split(':')[-1]))
