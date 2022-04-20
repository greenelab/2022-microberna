library(dplyr)
library(readr)

#species_string <- 's__Pseudomonas_aeruginosa'
species_string <- snakemake@wildcards[['gtdb_species']]
gtdb_species_string <- gsub("(s__.*?)_", "\\1 ", species_string)
gtdb_lineages <- read_csv(snakemake@input[['gtdb_lineages']])

gtdb_lineages_species <- gtdb_lineages %>%
  filter(species == gtdb_species_string)
write_csv(gtdb_lineages_species, snakemake@output[['accessions']])
