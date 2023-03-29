from promoter_seq import align_promoters, load_promoter_seq
import random
from multiprocessing import Pool
import pickle
from pyfaidx import Fasta

class BaselineCalculator(object):
    def __init__(self, prefix, export_prefix, fasta, seq_len, k=10000):
        self.fasta = Fasta(fasta, as_raw=False)
        self.canons = pickle.load(open('/data/db/import/save/'+export_prefix+'-canons.pkl', 'rb'))
        self.prefix = prefix
        self.seq_len = seq_len
        self.k = k

    def do(self):
        p=Pool(8)
        genes = [load_promoter_seq(r[0], self.fasta, self.seq_len).seq.upper() for r in self.canons.values()]
        self.scores = p.map(align_promoters, [(x,y) for x,y in zip(random.sample(genes, self.k), random.sample(genes, self.k))], chunksize=1000)
        pickle.dump(self.scores, open('/data/db/import/save/'+self.prefix+'-promoter-baseline-sim-%d.pkl'%self.seq_len, 'wb'))

baseline_hs_100 = BaselineCalculator('human', 'paralogs', r'/mnt/hddata/2/evgeny/TAU/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa', 100)
baseline_mm_100 = BaselineCalculator('mouse', 'mouse', r'/mnt/hddata/2/evgeny/TAU/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa', 100)
baseline_hs_300 = BaselineCalculator('human', 'paralogs', r'/mnt/hddata/2/evgeny/TAU/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa', 300)
baseline_mm_300 = BaselineCalculator('mouse', 'mouse', r'/mnt/hddata/2/evgeny/TAU/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa', 300)
baseline_hs_500 = BaselineCalculator('human', 'paralogs', r'/mnt/hddata/2/evgeny/TAU/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa', 500)
baseline_mm_500 = BaselineCalculator('mouse', 'mouse', r'/mnt/hddata/2/evgeny/TAU/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa', 500)
baseline_hs_1000 = BaselineCalculator('human', 'paralogs', r'/mnt/hddata/2/evgeny/TAU/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa', 1000)
baseline_mm_1000 = BaselineCalculator('mouse', 'mouse', r'/mnt/hddata/2/evgeny/TAU/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa', 1000)

if __name__ == '__main__':
    #baseline_hs_100.do()
    #baseline_mm_100.do()
    #baseline_hs_300.do()
    #baseline_mm_300.do()
    baseline_hs_500.do()
    baseline_mm_500.do()
    baseline_hs_1000.do()
    baseline_mm_1000.do()
