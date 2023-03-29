
def paralogs(db, fasta, specie):
    from overlap_calc import canons_by_gene
    from math import nan
    def f_or_n(f):
       try:
           return float(f)
       except ValueError:
           return nan

    import csv
    from collections import namedtuple
    species = set()
    branch_points = {}
    para = namedtuple('paralog', ['lgene', 'rgene', 'branch_point', 'species', 'dn', 'ds', 'sim'])
    i=0
    canons = canons_by_gene(db, fasta)
    with open('/data/db/within_species_paralogs_dnds_allcov.csv') as file:
        for p in (para(l['lgene'], l['rgene'], l['branch_point'], l['species'], f_or_n(l['dn']), f_or_n(l['ds']), f_or_n(l['sim'])) for l in csv.DictReader(file, fieldnames=['lgene', 'rgene', 'branch_point', 'species', 'dn', 'ds', 'sim'], quotechar='"')):
            if p.species != specie:
                continue
            i += 1
            if i % 1024 == 0:
                print(".",end='',flush=True)
            if not p.branch_point in branch_points:
                branch_points[p.branch_point] = []
            if ('gene:'+p.lgene) not in canons or ('gene:'+p.rgene) not in canons:
                continue
            lc = canons['gene:'+p.lgene]
            rc = canons['gene:'+p.rgene]
            l = (abs(lc[0].start-lc[0].end)+abs(rc[0].start-rc[0].end))/2
            dist = 300000000 if lc[0].chrom != rc[0].chrom or lc[0].strand != rc[0].strand else abs(lc[0].start-rc[0].start) if lc[0].strand == '+' else abs(lc[0].end-rc[0].end)
            branch_points[p.branch_point].append(('gene:'+p.lgene,'gene:'+p.rgene, p.dn, p.ds, l, dist, p.sim))
    return branch_points

def orthologs_by_gene(genes):
    import csv
    from collections import namedtuple
    species = set()
    o2os = { k:set() for k in genes }
    ortho = namedtuple('one2one', ['lgene', 'lspecies', 'rgene', 'rspecies'])
    with open('/data/db/one2one.csv') as file:
        for o in (ortho(l['gl.stable_id'], l['namel.name'], l['gr.stable_id'], l['namer.name']) for l in csv.DictReader(file, fieldnames=['gl.stable_id', 'namel.name', 'gr.stable_id', 'namer.name'], quotechar='"')):
            species.add(o.rspecies)
            if 'gene:'+o.lgene in o2os:
                o2os['gene:'+o.lgene].add(o.rspecies)
    return o2os

def ortholog_count_by_gene(genes):
    import csv
    from collections import namedtuple
    from pandas import DataFrame
    species = set()
    o2os = { k:{} for k in genes }
    ortho = namedtuple('one2one', ['lgene', 'lspecies', 'rgene', 'rspecies'])
    with open('/data/db/one2one.csv') as file:
        for o in (ortho(l['gl.stable_id'], l['namel.name'], l['gr.stable_id'], l['namer.name']) for l in csv.DictReader(file, fieldnames=['gl.stable_id', 'namel.name', 'gr.stable_id', 'namer.name'], quotechar='"')):
            species.add(o.rspecies)
            if 'gene:'+o.lgene in o2os:
                if o.rspecies in o2os['gene:'+o.lgene]:
                    o2os['gene:'+o.lgene][o.rspecies] += 1
                else:
                    o2os['gene:'+o.lgene][o.rspecies] = 1
    return DataFrame.from_records([(k,)+tuple((o2os[k][s]  if s in o2os[k] else 0 for s in species)) for k in o2os], columns=['gene']+list(species)).set_index('gene')
    

def write_paralog_sets(exp, fn):
    d = {}
    for i in exp.iloc:
        g1 = i.g1
        g2 = i.g2
        if g1 in d and g2 in d:
            d[g1]=d[g2]=d[g1].union(d[g2])
        elif g1 in d:
            d[g2]=d[g1]=d[g1].union(set([g2]))
        elif g2 in d:
            d[g1]=d[g2]=d[g2].union(set([g1]))
        else:
            d[g1]=d[g2]=set([g1, g2])
    with open(fn, 'w') as f:
        for i in d.keys():
            if d[i] != None:
                f.write(','.join(list(d[i]))+'\n')
                for j in d[i]:
                    d[j]=None
