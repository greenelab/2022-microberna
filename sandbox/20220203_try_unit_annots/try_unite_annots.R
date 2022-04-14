library(readr)
library(dplyr)
source("scripts/utils.R")

setwd("~/github/2022-microberna/")

# Explore whether annotations can be united across different information streams:
# 1. eggnog ortholog annotation of a pangenome
# 2. bakta annotation of a genome
# 3. RefSeq annotation of a genome
# 
# Test genome: GCF_010509575.1 is the accession number for Faecalibacterium prausnitzii A2-165,
# the strain that we used to analyze three benchmarking data sets.

# download test RefSeq annotation file

destfile <- "sandbox/20220203_try_unit_annots/GCF_010509575.1_ASM1050957v1_genomic.gff.gz"
url <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/010/509/575/GCF_010509575.1_ASM1050957v1/GCF_010509575.1_ASM1050957v1_genomic.gff.gz"
if (!file.exists(destfile)) {
  download.file(url, destfile, method="auto") 
}

# read in various annotation files that need to be united -----------------

refseq_gff <- read_gff("sandbox/20220203_try_unit_annots/GCF_010509575.1_ASM1050957v1_genomic.gff.gz")
bakta_gff <- read_gff("outputs/gtdb_genomes_bakta/s__Faecalibacterium_prausnitzii_C/GCF_010509575.1.gff3")
pangenome_clstrs <- read_cdhit_clstr("outputs/gtdb_genomes_annotated_comb/s__Faecalibacterium_prausnitzii_C_clustered_annotated_seqs.fa.clstr")

# label batka cols
bakta_gff <- bakta_gff %>%
  select(seqname, start, end, bakta_id = ID, bakta_name = Name)

# create a map from refseq features to bakta annotations
bakta_to_refseq <- refseq_gff %>%
  full_join(bakta_gff, by = c("seqname", "start", "end")) 

# create a map from bakta to pangenome clusters
bakta_to_pangenome <- bakta_gff %>%
  left_join(pangenome_clstrs, by = c("bakta_id" = "sequence_name")) %>%
  select(bakta_id, bakta_name, pangenome_representative = representative)

bakta_pangenome_refseq_map <- full_join(bakta_to_pangenome, bakta_to_refseq, by = c("bakta_id", "bakta_name"))

# tximport requires the following columns:
# TXNAME            GENEID
# pangenome_representative will always be the TXNAME, as that's what we quantified with salmon
# GENEID will vary depending on the desired output. 

# TODO: Add eggnog annotations to map; probably the most general use for the compendia