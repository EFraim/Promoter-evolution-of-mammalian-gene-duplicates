#!env python3
import sys
import re
gene_prefix = sys.argv[1]
for l in sys.stdin:
    m = re.findall(gene_prefix+'[0-9]+', l)
    if len(m) < 2:
        continue
    for x,y in zip(m[:-1], m[1:]):
        print('%s,%s'%(x,y))
