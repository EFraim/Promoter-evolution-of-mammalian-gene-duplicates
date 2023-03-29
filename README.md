# Promoter-evolution-of-mammalian-gene-duplicates
Scripts, figures and source code for "Promoter evolution of mammalian gene duplicates"

## Producing gene feature databases
1. Obtain [Homo_sapiens.GRCh38.98.gff3](https://ftp.ensembl.org/pub/release-98/gff3/homo_sapiens/Homo_sapiens.GRCh38.98.chr.gff3.gz)
2. Obtain [Mus_musculus.GRCm38.98.gff3](https://ftp.ensembl.org/pub/release-98/gff3/mus_musculus/Mus_musculus.GRCm38.98.chr.gff3.gz)
3.     import gffutils
       gffutils.create_db('Homo_sapiens.GRCh38.98.gff3', dbfn='homo.db', force=True, merge_strategy='merge')
       gffutils.create_db('Mus_musculus.GRCm38.98.gff3', dbfn='mouse.db', force=True, merge_strategy='merge')
       
## Calculating CpG overlaps 
