
def load_promoter_seq(g, fasta, len):
    if g.strand == '+':
        return fasta.get_seq(g.chrom, max(g.start-len, 0), max(g.start-len, 0)+min(len, g.start), rc=False)
    else:
        return fasta.get_seq(g.chrom, g.end, g.end+len, rc=True)

def load_polya_score(g, fasta, len):
    if g.strand == '+':
        seq=fasta.get_seq(g.chrom, max(g.end-len, g.start), g.end, rc=False).seq.upper()
    else:
        seq=fasta.get_seq(g.chrom, g.start, min(g.start+len, g.end), rc=True).seq.upper()
    cscore,maxscore=0,0
    for a in seq:
        cscore = cscore + 1 if a == 'A' else max(cscore-1, 0)
        maxscore = max(cscore, maxscore)
    return maxscore

def load_polya_kmer_score(g, fasta, len):
    if g.strand == '+':
        seq=fasta.get_seq(g.chrom, max(g.end-len, g.start), g.end, rc=False).seq.upper()
    else:
        seq=fasta.get_seq(g.chrom, g.start, min(g.start+len, g.end), rc=True).seq.upper()
    cscore=0
    for a, b, c, d, e in zip(seq[:-4], seq[1:-3], seq[2:-2], seq[3:-1], seq[4:]):
        if a == 'A' and b == 'A' and c == 'A' and d == 'A' and e == 'A':
            cscore = cscore + 1 
    return cscore

match_dic = {}
match_dic[('A', 'A')] = 2
match_dic[('A', 'C')] = -2
match_dic[('A', 'T')] = -2
match_dic[('A', 'G')] = -1
match_dic[('C', 'A')] = -2
match_dic[('C', 'C')] = 2
match_dic[('C', 'T')] = -1
match_dic[('C', 'G')] = -2
match_dic[('G', 'A')] = -1
match_dic[('G', 'C')] = -2
match_dic[('G', 'T')] = -2
match_dic[('G', 'G')] = 2
match_dic[('T', 'A')] = -2
match_dic[('T', 'C')] = -1
match_dic[('T', 'T')] = 2
match_dic[('T', 'G')] = -2
match_dic[('N', 'A')] = -2
match_dic[('N', 'C')] = -2
match_dic[('N', 'T')] = -2
match_dic[('N', 'G')] = -2
match_dic[('A', 'N')] = -2
match_dic[('C', 'N')] = -2
match_dic[('T', 'N')] = -2
match_dic[('G', 'N')] = -2
match_dic[('N', 'N')] = 0
    
#cnt=0
def align_promoters(seq_g):
    from Bio import pairwise2
#    global cnt
#    cnt += 1
#    #if cnt % 16384 == 0:
#    print('.', end='', flush=True)
#    print(seq_g[0])
#    print(seq_g[1])
    return pairwise2.align.localds(seq_g[0], seq_g[1], match_dic, -8, -0.01)[0]

import pickle
from pyfaidx import Fasta
import Bio.Align
import numpy as np
from multiprocessing import Pool
import multiprocessing
import logging
multiprocessing.log_to_stderr().setLevel(logging.DEBUG)

class PromoterSeqAnalysis(object):
    def __init__(self, prefix, export_prefix, canons_prefix, fasta, seq_len, polya_len):
        self.fasta = Fasta(fasta, as_raw=False)
        self.canons = pickle.load(open('/data/db/import/save/'+canons_prefix+'-canons.pkl', 'rb'))
        self.export = pickle.load(open('/data/db/import/save/'+export_prefix+'-export.pkl', 'rb'))
        self.prefix = prefix
        self.seq_len = seq_len
        self.polya_len = polya_len

    def do(self):
        print('Processing %s(%d)'%(self.prefix, self.seq_len))
        p=Pool(8)
        genes = [(load_promoter_seq(self.canons['gene:'+r['g1']][0], self.fasta, self.seq_len).seq.upper(),
                  load_promoter_seq(self.canons['gene:'+r['g2']][0], self.fasta, self.seq_len).seq.upper()) for r in self.export.iloc]
        polya = [(load_polya_score(self.canons['gene:'+r['g1']][0], self.fasta, self.polya_len),
                  load_polya_score(self.canons['gene:'+r['g2']][0], self.fasta, self.polya_len)) for r in self.export.iloc]
        polyak = [(load_polya_kmer_score(self.canons['gene:'+r['g1']][0], self.fasta, self.polya_len),
                  load_polya_kmer_score(self.canons['gene:'+r['g2']][0], self.fasta, self.polya_len)) for r in self.export.iloc]
        alignments = p.map(align_promoters, genes, chunksize=16384)
        #alignments = [align_promoters(g) for g in genes[38*16384+14914:]] #p.map(align_promoters, genes, chunksize=16384)
        self.alignments = list(zip(((r['g1'], r['g2']) for r in self.export.iloc), alignments))
        self.export = self.export.assign(**{'PromPresScore%s'%self.seq_len: [a[2] for a in alignments],
                                          'PromPresTssMinus%s'%self.seq_len: [min(len(a[0][a[4]:])-a[0][a[4]:].count('-'), len(a[1][a[4]:])-a[1][a[4]:].count('-')) for a in alignments],
                                          'PromPresShift%s'%self.seq_len: [abs(len(a[0][a[4]:])-a[0][a[4]:].count('-')-len(a[1][a[4]:])-a[1][a[4]:].count('-')) for a in alignments],
                                          'PolyA1%s'%self.seq_len: [p[0] for p in polya], 'PolyA2%s'%self.seq_len: [p[1] for p in polya],
                                          'PolyAK1%s'%self.seq_len: [p[0] for p in polyak], 'PolyAK2%s'%self.seq_len: [p[1] for p in polyak]})
        self.export.assign = self.export.assign(trans1 = [v[0].id.split(':')[-1] for v in ('gene:'+self.export['g1']).map(self.canons)],
                                                trans2 = [v[0].id.split(':')[-1] for v in ('gene:'+self.export['g2']).map(self.canons)])
        pickle.dump(self.export, open('/data/db/import/save/'+self.prefix+'-export.pkl', 'wb'))

class PromoterOrthoSeqAnalysis(object):
    def __init__(self, lspecies, rspecies, lfasta, rfasta, ldb, rdb, prefix, seq_len):
        from gffutils import FeatureDB
        from csv import DictReader
        self.lfasta = Fasta(lfasta, as_raw=False)
        self.rfasta = Fasta(rfasta, as_raw=False)
        with open('/data/db/one2one.csv') as file:
            self.orthos = [(l['gl.stable_id'], l['gr.stable_id']) for l in DictReader(file, fieldnames=['gl.stable_id', 'namel.name', 'gr.stable_id', 'namer.name'], quotechar='\"') if l['namel.name'] == lspecies and l['namer.name'] == rspecies]
        self.ldb = FeatureDB(ldb)
        self.rdb = FeatureDB(rdb)
        self.prefix = prefix
        self.seq_len = seq_len

    def do(self):
        from pandas import DataFrame
        p=Pool(8)
        genes = [(load_promoter_seq(self.ldb['gene:'+r[0]], self.lfasta, self.seq_len).seq.upper(),
                  load_promoter_seq(self.rdb['gene:'+r[1]], self.rfasta, self.seq_len).seq.upper()) for r in self.orthos]
        alignments = p.map(align_promoters, genes, chunksize=16384)
        self.alignments = list(zip(self.orthos, alignments))
        self.export = DataFrame.from_records({'g1': a[0][0], 'g2': a[0][1], 'PromPresScore%s'%self.seq_len: a[1][2],
                                'PromPresTssMinus%s'%self.seq_len: min(len(a[1][0][a[1][4]:])-a[1][0][a[1][4]:].count('-'), len(a[1][1][a[1][4]:])-a[1][1][a[1][4]:].count('-')),
                                'PromPresShift%s'%self.seq_len: abs(len(a[1][0][a[1][4]:])-a[1][0][a[1][4]:].count('-')-len(a[1][1][a[1][4]:])-a[1][1][a[1][4]:].count('-'))} for a in self.alignments)
        pickle.dump(self.export, open('/data/db/import/save/'+self.prefix+'-ortho-sim-local.pkl', 'wb'))
        self.export.to_csv('/data/db/import/save/'+self.prefix+'-ortho-sim-local.csv')

        
# aligner_hs100 = PromoterSeqAnalysis('human-sim-local', 'paralogs', 'paralogs', r'/mnt/hddata/2/evgeny/TAU/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa', 100, 70)
# aligner_mm100 = PromoterSeqAnalysis('mouse-sim-local', 'mouse', 'mouse', r'/mnt/hddata/2/evgeny/TAU/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa', 100, 70)
# if __name__ == '__main__':
#     aligner_hs100.do()
#     aligner_mm100.do()
# aligner_hs300 = PromoterSeqAnalysis('human-sim-local', 'human-sim-local', 'paralogs', r'/mnt/hddata/2/evgeny/TAU/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa', 300, 70)
# aligner_mm300 = PromoterSeqAnalysis('mouse-sim-local', 'mouse-sim-local', 'mouse', r'/mnt/hddata/2/evgeny/TAU/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa', 300, 70)
# if __name__ == '__main__':
#     aligner_hs300.do()
#     aligner_mm300.do()
#aligner_hs500 = PromoterSeqAnalysis('human-sim-local', 'human-sim-local', 'paralogs', r'/mnt/hddata/2/evgeny/TAU/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa', 500, 70)
aligner_mm500 = PromoterSeqAnalysis('mouse-sim-local', 'mouse-sim-local', 'mouse', r'/mnt/hddata/2/evgeny/TAU/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa', 500, 70)
if __name__ == '__main__':
    #aligner_hs500.do()
    aligner_mm500.do()
aligner_hs1000 = PromoterSeqAnalysis('human-sim-local', 'human-sim-local', 'paralogs', r'/mnt/hddata/2/evgeny/TAU/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa', 1000, 70)
aligner_mm1000 = PromoterSeqAnalysis('mouse-sim-local', 'mouse-sim-local', 'mouse', r'/mnt/hddata/2/evgeny/TAU/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa', 1000, 70)
#aligner_hm = PromoterOrthoSeqAnalysis('Homo sapiens', 'Mus musculus', r'/mnt/hddata/2/evgeny/TAU/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa', r'/mnt/hddata/2/evgeny/TAU/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa', 'homo.db', 'mouse.db', 'human-mouse', 300)
if __name__ == '__main__':
    aligner_hs1000.do()
    aligner_mm1000.do()
    #aligner_hm.do()
