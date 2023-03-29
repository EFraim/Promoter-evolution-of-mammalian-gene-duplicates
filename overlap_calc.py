
def canonical_transcripts(db, fasta_filename):
    import pyfaidx
    fasta = pyfaidx.Fasta(fasta_filename, as_raw=False)
    best_transcript = None
    for gene in db.features_of_type('gene'):
        # exons_list will contain (CDS_length, total_length, transcript, [exons]) tuples.
        exon_list = []
        #print(gene)
        for ti, transcript in enumerate(db.children(gene, level=1)):
            if transcript.attributes['biotype'][0] != 'protein_coding':
                continue
            total_len = 0
            exons = list(db.children(transcript, level=1))
            for exon in exons:
                exon_length = len(exon)
                if exon.featuretype == 'exon':
                    total_len += exon_length

            exon_list.append((total_len, transcript, [e for e in exons if e.featuretype in ['exon']] ))

        # Pick the longest transcript among protein_coding ones with the lowest tls
        #print(exon_list)
        if len(exon_list) < 1:
            continue
        best = sorted(exon_list, key=lambda x: (x[1].attributes['transcript_support_level'][0].split(' ')[0] if 'transcript_support_level' in x[1].attributes and len(x[1].attributes['transcript_support_level'])>0 else 'NA', -x[0]))[0]

        canonical_exons = best[-1]
        transcript = best[-2]
        seqs = [i.sequence(fasta) for i in sorted(canonical_exons, key=lambda x: x.start, reverse=transcript.strand != '+')]
        first_exons=sorted(canonical_exons, key=lambda x: x.start, reverse=transcript.strand != '+')[0:2] if len(canonical_exons)>=2 else None
        yield transcript, ''.join(seqs), (first_exons[1].start-first_exons[0].end if transcript.strand == '+' else first_exons[0].start - first_exons[1].end) if len(canonical_exons) >= 2 else 0, len(canonical_exons)-1

def canonical_exons(db, fasta_filename):
    import pyfaidx
    fasta = pyfaidx.Fasta(fasta_filename, as_raw=False)
    best_transcript = None
    for gene in db.features_of_type('gene'):
        # exons_list will contain (CDS_length, total_length, transcript, [exons]) tuples.
        exon_list = []
        #print(gene)
        for ti, transcript in enumerate(db.children(gene, level=1)):
            if transcript.attributes['biotype'][0] != 'protein_coding':
                continue
            total_len = 0
            exons = list(db.children(transcript, level=1))
            for exon in exons:
                exon_length = len(exon)
                if exon.featuretype == 'exon':
                    total_len += exon_length

            exon_list.append((total_len, transcript, [e for e in exons if e.featuretype in ['exon']] ))

        # Pick the longest transcript among protein_coding ones with the lowest tls
        #print(exon_list)
        if len(exon_list) < 1:
            continue
        best = sorted(exon_list, key=lambda x: (x[1].attributes['transcript_support_level'][0].split(' ')[0] if 'transcript_support_level' in x[1].attributes and len(x[1].attributes['transcript_support_level'])>0 else 'NA', -x[0]))[0]

        canonical_exons = sorted(best[-1], key=lambda x: x.start)
        transcript = best[-2]
        seqs = [i.sequence(fasta) for i in sorted(canonical_exons, key=lambda x: x.start, reverse=transcript.strand != '+')]
        yield gene.id, ''.join(seqs), canonical_exons

def canonical_coding_exons(db, fasta_filename):
    import pyfaidx
    fasta = pyfaidx.Fasta(fasta_filename, as_raw=False)
    best_transcript = None
    for gene in db.features_of_type('gene'):
        # exons_list will contain (CDS_length, total_length, transcript, [exons]) tuples.
        exon_list = []
        #print(gene)
        for ti, transcript in enumerate(db.children(gene, level=1)):
            if transcript.attributes['biotype'][0] != 'protein_coding':
                continue
            total_len = 0
            exons = list(db.children(transcript, level=1))
            for exon in exons:
                exon_length = len(exon)
                if exon.featuretype == 'CDS':
                    total_len += exon_length

            exon_list.append((total_len, transcript, [e for e in exons if e.featuretype in ['CDS']] ))

        # Pick the longest transcript among protein_coding ones with the lowest tls
        #print(exon_list)
        if len(exon_list) < 1:
            continue
        best = sorted(exon_list, key=lambda x: (x[1].attributes['transcript_support_level'][0].split(' ')[0] if 'transcript_support_level' in x[1].attributes and len(x[1].attributes['transcript_support_level'])>0 else 'NA', -x[0]))[0]

        canonical_exons = sorted(best[-1], key=lambda x: x.start)
        transcript = best[-2]
        seqs = [i.sequence(fasta) for i in sorted(canonical_exons, key=lambda x: x.start, reverse=transcript.strand != '+')]
        yield gene.id, ''.join(seqs), canonical_exons

def canonical_transcripts_nc(db, fasta_filename):
    import pyfaidx
    fasta = pyfaidx.Fasta(fasta_filename, as_raw=False)
    for gene in db.features_of_type(['gene', 'ncRNA_gene', 'snoRNA_gene', 'gRNA_gene','RNase_P_RNA_gene','tmRNA_gene','SRP_RNA_gene','tRNA_gene','enzymatic_RNA_gene','lncRNA_gene','rRNA_gene','RNase_MRP_RNA_gene','snRNA_gene','telomerase_RNA_gene','miRNA_gene','piRNA_gene','scRNA_gene', 'pseudogene']):

        # exons_list will contain (CDS_length, total_length, transcript, [exons]) tuples.
        exon_list = []
        for ti, transcript in enumerate(db.children(gene, level=1)):
            cds_len = 0
            total_len = 0
            exons = list(db.children(transcript, level=1))
            for exon in exons:
                exon_length = len(exon)
                if exon.featuretype == 'CDS':
                    cds_len += exon_length
                total_len += exon_length

            exon_list.append((cds_len, total_len, transcript, exons if cds_len == 0 else [e for e in exons if e.featuretype in ['exon']] ))

        # If we have CDS, then use the longest coding transcript
        if max(i[0] for i in exon_list) > 0:
            best = sorted(exon_list, key=lambda x: x[0], reverse=True)[0]
        # Otherwise, just choose the longest
        else:
            best = sorted(exon_list, key=lambda x: x[1])[0]

        #print(best)

        canonical_exons = best[-1]
        transcript = best[-2]
        seqs = [i.sequence(fasta) for i in sorted(canonical_exons, key=lambda x: x.start, reverse=transcript.strand != '+')]
        first_exons=sorted(canonical_exons, key=lambda x: x.start, reverse=transcript.strand != '+')[0:2] if len(canonical_exons)>=2 else None
        yield transcript, ''.join(seqs), (first_exons[1].start-first_exons[0].end if transcript.strand == '+' else first_exons[0].start - first_exons[1].end) if len(canonical_exons) >= 2 else 0, len(canonical_exons)-1


def canons_by_gene(db, fa, nc=False):
    import gffutils
    canons = list(g for g in (canonical_transcripts_nc(db, fa) if nc else canonical_transcripts(db, fa)))
    if not nc:
        canons = [g for g in canons if g[0].attributes['biotype'][0] == 'protein_coding']
    return dict(((g[0].attributes['Parent'][0], g) for g in canons))

def load_cpg_gff(db):
    return [i for i in db.features_of_type('biological_region') if 'cpg' in i.attributes['logic_name']]

def load_cpg_tsv(fname):
    import csv
    from collections import namedtuple
    CPG = namedtuple('CPG', ['chrom', 'start', 'end'])
    with open(fname) as file:
        return [CPG(l['chrom'][3:].upper(), int(l['chromStart']), int(l['chromEnd'])) for l in csv.DictReader(file, delimiter='\t')]

def load_cpg_bed(fname):
    import csv
    from collections import namedtuple
    CPG = namedtuple('CPG', ['chrom', 'start', 'end'])
    with open(fname) as file:
        return [CPG(l['chrom'][3:].upper(), int(l['chromStart']), int(l['chromEnd'])) for l in csv.DictReader(file, delimiter='\t', fieldnames=['chrom', 'chromStart', 'chromEnd'])]
    

def first_intron(db):
    for gene in db.features_of_type('gene'):

        # exons_list will contain (CDS_length, total_length, transcript, [exons]) tuples.
        exon_list = []
        for ti, transcript in enumerate(db.children(gene, level=1)):
            cds_len = 0
            total_len = 0
            exons = list(db.children(transcript, level=1))
            for exon in exons:
                exon_length = len(exon)
                if exon.featuretype == 'CDS':
                    cds_len += exon_length
                total_len += exon_length

            exon_list.append((cds_len, total_len, transcript, exons if cds_len == 0 else [e for e in exons if e.featuretype in ['CDS', 'five_prime_UTR', 'three_prime_UTR']] ))

        # If we have CDS, then use the longest coding transcript
        if max(i[0] for i in exon_list) > 0:
            best = sorted(exon_list, key=lambda x: x[0], reverse=True)[0]
        # Otherwise, just choose the longest
        else:
            best = sorted(exon_list, key=lambda x: x[1])[0]

        #print(best)

        canonical_exons = best[-1]
        if len(canonical_exons) < 2:
            yield gene, 0
        else:
            first_exons=sorted(canonical_exons, key=lambda x: x.start, reverse=transcript.strand != '+')[0:2]
            yield gene, first_exons[1].start-first_exons[0].end if transcript.strand == '+' else first_exons[0].start - first_exons[1].end


def cpg_overlap(canons, cpgs, pre_tss=100, post_tss=300):
    gene_overlap = {}
    for gene in canons:
        promoter = [0]*(pre_tss+post_tss+1)
        strandpos = gene[0].strand == '+'
        for island in cpgs:
            if island.chrom != gene[0].chrom:
                continue
            tss = gene[0].start if strandpos else gene[0].end
            for i in range(max(island.start, tss-(pre_tss if strandpos else post_tss))-tss+(pre_tss if strandpos else post_tss),
                           min(island.end, tss+(post_tss if strandpos else pre_tss)+1)-tss+(pre_tss if strandpos else post_tss)):
                promoter[i] = 1
        gene_overlap[gene[0].attributes['Parent'][0]] = sum(promoter)/len(promoter)
    return gene_overlap

def split_overlaps(gene_overlap):
    return [g for g in gene_overlap.keys() if gene_overlap[g] < 0.5], [g for g in gene_overlap.keys() if gene_overlap[g] > 0.5]

def canon_ints(canons, pre_tss=100, post_tss=300):
    from pybedtools import cbedtools
    ints = []
    for gene in canons:
        strandpos = gene[0].strand == '+'
        tss = gene[0].start if strandpos else gene[0].end
        ints.append(cbedtools.Interval(gene[0].chrom if gene[0].chrom[:3]=='chr' else 'chr'+gene[0].chrom,
                                       max(tss-(pre_tss if strandpos else post_tss)-1, 0),
                                       tss+(post_tss if strandpos else pre_tss)+1-1,
                                       name=gene[0].attributes['Parent'][0], strand=gene[0].strand))
    return ints


def ordered_cpgs(cpgs):
    from sortedcontainers import SortedDict
    res = {}
    for chrom in set([i.chrom for i in cpgs]):
        res[chrom] = SortedDict([(i.start, i.end) for i in cpgs if i.chrom == chrom])
    return res

def nmi_tss(chrom, tss, d):
  if not chrom in d:
    return False
  bp = d[chrom].bisect_left(tss)
  if bp < len(d[chrom]) and tss >= d[chrom].peekitem(bp)[0]-1000 and tss <= d[chrom].peekitem(bp)[1]+1000:
    return True
  if bp > 0 and tss >= d[chrom].peekitem(bp-1)[0]-1000 and tss <= d[chrom].peekitem(bp-1)[1]+1000:
    return True
  return False
