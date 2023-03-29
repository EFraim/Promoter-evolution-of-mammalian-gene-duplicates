from overlap_calc import *

def iter_cigar(cigar):
    i=0
    while i<len(cigar):
        if cigar[i] == 'D' or cigar[i] == 'M':
            yield cigar[i]
            i += 1
        else:
            j=i
            while cigar[j].isdigit():
                j += 1
            for k in range(0, int(cigar[i:j])):
                yield cigar[j]
            i = j+1
                
def adjusted_exons(exon, cigs):
    from itertools import accumulate
    pts = list(accumulate([e.end-e.start+1 for e in (exon if exon[0].strand == '+' else exon[::-1])][:-1]))
    offset = 0
    consumed = 0
    idx=0
    for c in iter_cigar(cigs):
        while idx < len(pts) and offset + consumed >= pts[idx]:
            pts[idx] += offset
            idx += 1
        if c == 'M':
            consumed += 3
        elif c == 'D':
            offset += 3
    #print(pts)
    return pts
            
def consistent_exons(exo1, exo2, dist):
    return 1.0 if len(exo1)==0 and len(exo2)==0 else 0.0 if len(exo1)==0 or len(exo2)==0 else min(
        sum(min(abs(x-y) for y in exo2) <= dist for x in exo1)/len(exo1),
        sum(min(abs(x-y) for y in exo1) <= dist for x in exo2)/len(exo2))

if __name__ == "__main__":
    import gffutils
    import pickle
    db = gffutils.FeatureDB('mouse.db')
    exons={ d[0][5:]:d[2] for d in canonical_exons(db, r'/mnt/hddata/2/evgeny/TAU/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa')}
    coding_exons={ d[0][5:]:d[2] for d in canonical_coding_exons(db, r'/mnt/hddata/2/evgeny/TAU/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa')}
    cigs = {}
    from csv import reader
    for r in reader(open('/data/db/import/mouse_paralog_cigarlines.csv', 'r')):
        if not r[0] in cigs:
            cigs[r[0]] = ((r[7], r[3]),)
        else:
            cigs[r[0]] = ((r[7], r[3]),cigs[r[0]][0])

    cigs = {(v[0][0], v[1][0]) : (v[0][1], v[1][1]) for v in cigs.values()}
    cigs.update({(k[1], k[0]) : (v[1],v[0]) for (k,v) in cigs.items()})

    export=pickle.load(open('/data/db/import/save/mouse-sim-local-export.pkl', 'rb'))

    exonsg1 = [exons[r['g1']] for r in export.iloc]
    exonsg2 = [exons[r['g2']] for r in export.iloc]
    conexons=[consistent_exons(adjusted_exons(coding_exons[r['g1']], cigs[(r['g1'], r['g2'])][0]),
                               adjusted_exons(coding_exons[r['g2']], cigs[(r['g1'], r['g2'])][1]), 10) for r in export.iloc]
    export = export.assign(conexons=conexons)
    export = export.assign(exonsg1=[len(x) for x in exonsg1])
    export = export.assign(exonsg2=[len(x) for x in exonsg2])
    export = export.assign(retro = [ -1 if len(r[0])>=3 and len(r[1])==1 else 1 if len(r[1])>=3 and len(r[0])==1 else 0 for r in  zip(exonsg1, exonsg2)])
    export = export.assign(segmental=[r[0] >= 0.8 and len(r[1]) >= 3 and len(r[2]) >= 3 and len(r[1])/len(r[2]) >= 0.8 and len(r[2])/len(r[1]) >= 0.8 for r in zip(conexons, exonsg1, exonsg2)])
    export[export['g1']<export['g2']].to_csv('mouse-paralogs-exocon.csv')
