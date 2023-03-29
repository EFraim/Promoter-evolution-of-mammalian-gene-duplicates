
from pandas import *
from overlap_calc import *
import gffutils
import numpy as np
import pickle
import re
import matplotlib.pyplot as plt
from csv import reader, DictReader
CHROMLEN=200000000



db = gffutils.FeatureDB('mouse.db')
mouse_canons = list(g for g in gffutils.helpers.canonical_transcripts(db, r'/mnt/hddata/2/evgeny/TAU/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa') if g[0].attributes['biotype'][0] == 'protein_coding')
cpgs_mouse = load_cpg_gff(db)
overlap = cpg_overlap(mouse_canons, cpgs_mouse, pre_tss=300, post_tss=100)
cpg_less, cpg_full = [set(s) for s in split_overlaps(overlap)]
canons = canons_by_gene(db, r'/mnt/hddata/2/evgeny/TAU/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa')
facthist = { c : np.zeros(CHROMLEN, dtype='uint16') for c in set(canon[0].chrom for canon in canons.values()) }

DCids = []
for prom in DictReader(open('/data/db/import/mouse_factor_full_QC.txt', 'r'), delimiter='\t'):
    if re.fullmatch("CTCF|RAD21|REST|EP300|.*POLR.*", prom['Factor']):
        continue
    for inst in reader(open('/data/db/import/mouse_factor/'+prom['DCid']+'_sort_peaks.narrowPeak.bed', 'r'), delimiter='\t'):
        if not inst[0].startswith('chr') or inst[0][3:] not in facthist: #It is not on chromosomes with coding genes
            continue
        facthist[inst[0][3:]][int(inst[1]):int(inst[2])+1] += 1

PRE=1000
POST=1000

cpglesshist = np.zeros(PRE+POST, dtype='uint32')
cpgfullhist = np.zeros(PRE+POST, dtype='uint32')
cpglesserr = np.zeros(PRE+POST, dtype='uint32')
cpgfullerr = np.zeros(PRE+POST, dtype='uint32')
#for 
for cset, h in ((cpg_less, cpglesshist), (cpg_full, cpgfullhist)):
    for ckey in cset:
        c = canons[ckey][0]
        if c.strand == '+':
            h[max(0, PRE-c.start):min(POST+PRE, CHROMLEN-max(0, c.start-PRE))] += facthist[c.chrom][max(0, c.start-PRE):min(max(0, c.start-PRE)+POST+PRE-max(0, PRE-c.start), CHROMLEN)]
        else:
            h[max(0, POST-c.end):min(POST+PRE, CHROMLEN-max(0, c.end-POST))] += facthist[c.chrom][max(0, c.end-POST):min(max(0, c.end-POST)+POST+PRE-max(0, POST-c.end), CHROMLEN)][::-1]

cpglessavg = cpglesshist/float(len(cpg_less))
cpgfullavg = cpgfullhist/float(len(cpg_full))
cpglesserr = np.zeros(PRE+POST, dtype='float')
cpgfullerr = np.zeros(PRE+POST, dtype='float')

for cset, h, mean in ((cpg_less, cpglesserr, cpglessavg), (cpg_full, cpgfullerr, cpgfullavg)):
    for ckey in cset:
        c = canons[ckey][0]
        if c.strand == '+':
            h[max(0, PRE-c.start):min(POST+PRE, CHROMLEN-max(0, c.start-PRE))] += np.power(facthist[c.chrom][max(0, c.start-PRE):min(max(0, c.start-PRE)+POST+PRE-max(0, PRE-c.start), CHROMLEN)]-mean, 2.0)
        else:
            h[max(0, POST-c.end):min(POST+PRE, CHROMLEN-max(0, c.end-POST))] += np.power(facthist[c.chrom][max(0, c.end-POST):min(max(0, c.end-POST)+POST+PRE-max(0, POST-c.end), CHROMLEN)][::-1]-mean, 2.0)

cpglesserr = np.power(cpglesserr/len(cpg_less), 0.5)
cpgfullerr = np.power(cpgfullerr/len(cpg_full), 0.5)

plt.errorbar(range(-PRE, POST), cpglessavg, yerr=cpglesserr, errorevery=50, label='CpG-less', capsize=5)
plt.errorbar(range(-PRE, POST), cpgfullavg, yerr=cpgfullerr, errorevery=50, label='CpG-full', capsize=5)
plt.legend()
plt.savefig('mouse-promoter-err.png')

