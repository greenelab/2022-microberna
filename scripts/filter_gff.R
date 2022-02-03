library(rtracklayer)

# import gff annotation files and filter to remove contig regions, then rewrite in GFF format
gff <- readGFF(snakemake@input[['gff']], filter = list(type=c("CDS", "ncRNA", "tRNA", "rRNA", "tmRNA")))
export(gff, snakemake@output[['gff']])
