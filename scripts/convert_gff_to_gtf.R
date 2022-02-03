library(rtracklayer)

annot <- import(snakemake@input[['gff']])
export(annot, snakemake@output[['gtf']], "gtf")
