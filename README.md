# Promoter-evolution-of-mammalian-gene-duplicates
Scripts, figures and source code for "Promoter evolution of mammalian gene duplicates"

## Producing gene feature databases
1. Obtain [Homo_sapiens.GRCh38.98.gff3](https://ftp.ensembl.org/pub/release-98/gff3/homo_sapiens/Homo_sapiens.GRCh38.98.chr.gff3.gz)
2. Obtain [Mus_musculus.GRCm38.98.gff3](https://ftp.ensembl.org/pub/release-98/gff3/mus_musculus/Mus_musculus.GRCm38.98.chr.gff3.gz)
3.     import gffutils
       gffutils.create_db('Homo_sapiens.GRCh38.98.gff3', dbfn='homo.db', force=True, merge_strategy='merge')
       gffutils.create_db('Mus_musculus.GRCm38.98.gff3', dbfn='mouse.db', force=True, merge_strategy='merge')
       
## Calculating CpG overlaps 
This is taken care of by `calc_overlaps.py` script. It needs [Homo sapiens assembly](https://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz) and [Mus musculus assembly](https://ftp.ensembl.org/pub/release-98/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz)

## Paralogue analysis step
`Paralogs.ipynb` and `Mouse-Ensembl.ipynb` need to be evaluated.

### Outputs of this step
* `paralogs-export.pkl`
* `mouse-export.pkl`
* `human-genes.csv`
* `mouse-genes.csv`

## Calculating promoter sequence similarity
PSS in our paper is based on local sequence alignment, as described in Methods. It can be computed by using the `promoter_seq.py` script and the `promoter_baselines.py` script to compute the baseline for randomly selected pairs.

### Inputs and outputs of this step
* `paralogs-export.pkl` -> `human-sim-local-export.pkl`
* `mouse-export.pkl` -> `mouse-sim-local-export.pkl`

## Calculating exon congruence between paralogs
This step requires cigarlines for the paralogs. Those can be obtained from Ensembl database (release 98) but to save the trouble the precomputed version is available here in files `human_paralog_cigarlines.csv.xz` and `mouse_paralog_cigarlines.csv.xz`.

Those two files need to be extracted and then the `exons_alignment-human.py` and `exons_alignment-mouse.py` need to be executed.

### Inputs and outputs of this step
* `human-sim-local-export.pkl` -> `paralogs-export-exocon.pkl`, `human-paralogs-exocon.csv`, `human-paralogs-exocon-utr.csv`
* `mouse-sim-local-export.pkl` -> `mouse-paralogs-exocon.csv`

Note: the distinction between `human-paralogs-exocon.csv` and `human-paralogs-exocon-utr.csv` is due to some annotations having no UTR in this version of Ensembl for homo sapiens. Those turned out to be badly annotated genes, which had to be filtered not to introduce unacceptable noise into further analysis.

## Orthologs analysis
This is done by evaluating `Orthologs.ipynb`

## Promoter sequence similarity/mechanism analysis
This is done by evaluating `Promoter-Seq-Similarity.ipynb` and `Duplication-Mechanism.ipynb`

