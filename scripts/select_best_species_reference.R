library(readr)
library(dplyr)

#gtdb_lineages <- read_csv("inputs/gtdb-rs202.taxonomy.v2.csv")
gtdb_lineages <- read_csv(snakemake@input[['gtdb_lineages']])

#gather <- read_csv("outputs/rnaseq_sourmash_gather/ERX3558803_gtdb_k31.csv") %>%
gather <- read_csv(snakemake@input[['gather']]) %>%
  filter(gather_result_rank == 0) %>%
  select(name, query_name) %>%
  mutate(genome_accession = gsub(" .*", "", name)) %>%
  left_join(gtdb_lineages, by = c("genome_accession" = "ident")) %>%
  mutate(species_no_space = gsub(' ', '_', species)) %>%
  mutate(sra_to_ref_species = paste0(species_no_space, "-",  query_name))

write_csv(gather, snakemake@output[['sra_to_ref_species']])

