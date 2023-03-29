import gffutils
import pickle
from overlap_calc import *

db = gffutils.FeatureDB('homo.db')
mouse_db = gffutils.FeatureDB('mouse.db')
canons = list(g for g in canonical_transcripts(db, r'/mnt/hddata/2/evgeny/TAU/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa'))
mouse_canons = list(g for g in canonical_transcripts(mouse_db, r'/mnt/hddata/2/evgeny/TAU/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa'))

cpgs_compara = load_cpg_gff(db)
cpgs_ucsc = load_cpg_tsv("./cpg_islands_ucsc.tsv")
cpgs_mouse = load_cpg_gff(mouse_db)
cpgs_mouse_liver = load_cpg_bed('/home/evgeny/tau/gen_cpg/mm_liver_nmi.bed')
hs_ordered = ordered_cpgs(cpgs_ucsc)
mm_ordered = ordered_cpgs(cpgs_mouse_liver)
hs_cpgfull_elife={c[0].attributes['Parent'][0] : 1.0 if nmi_tss(c[0].chrom, c[0].start if c[0].strand=='+' else c[0].end, hs_ordered) else 0.0 for c in canons}
mm_cpgfull_elife={c[0].attributes['Parent'][0] : 1.0 if nmi_tss(c[0].chrom, c[0].start if c[0].strand=='+' else c[0].end, mm_ordered) else 0.0 for c in mouse_canons}

for pre in (300, 500, 1000, 2000):
    compara_overlap = cpg_overlap(canons, cpgs_compara, pre_tss=pre, post_tss=100)
    ucsc_overlap = cpg_overlap(canons, cpgs_ucsc, pre_tss=pre, post_tss=100)
    mouse_overlap = cpg_overlap(mouse_canons, cpgs_mouse, pre_tss=pre, post_tss=100)
    mouse_overlap_liver = cpg_overlap(mouse_canons, cpgs_mouse_liver, pre_tss=pre, post_tss=100)
    pickle.dump(compara_overlap, open('compara_overlaps_%d.pkl'%(pre,), 'wb'))
    pickle.dump(ucsc_overlap, open('ucsc_overlaps_%d.pkl'%(pre,), 'wb'))
    pickle.dump(mouse_overlap, open('mouse_overlaps_ensembl_%d.pkl'%(pre,), 'wb'))
    pickle.dump(mouse_overlap_liver, open('mouse_overlaps_liver_%d.pkl'%(pre,), 'wb'))

pickle.dump(hs_cpgfull_elife, open('hs_elife_overlaps.pkl', 'wb'))
pickle.dump(mm_cpgfull_elife, open('mm_elife_overlaps.pkl', 'wb'))
