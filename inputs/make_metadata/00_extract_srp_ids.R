library(dplyr)
library(readr)
setwd("~/github/2022-microberna/")

# This script reads in the SRA run table that annotates our RNAseq data of interest,
# filters that table to remove sequencing data sets we don't want to process
# (ex. chip seq, size fractionation selection), and writes the sample accessions
# to a file to be used by the next script to get runinfo tables from the ENA.

# annotate JGI accessions that will fail -----------------------------------

# when I navigate to the accessions that fail download, I see the banner
# "No public data has been made available in this project yet. Awaiting submission and/or validation of data."
# Accessions that failed:
jgi_sras_that_fail <- c("SRP235317", "SRP234344", "SRP234343", "SRP234342", 
                        "SRP234341", "SRP234340", "SRP234339", 'SRP234338', 
                        'SRP234337', 'SRP234336', 'SRP234335', 'SRP234334', 
                        'SRP234333', 'SRP234332', 'SRP234331', 'SRP234330',
                        'SRP234329', 'SRP234328', 'SRP234327', 'SRP234326', 
                        'SRP234325', 'SRP234324', 'SRP234323', 'SRP234322', 
                        'SRP234321', 'SRP234320', 'SRP234319', 'SRP234318',
                        'SRP231538', 'SRP231537', 'SRP231536', 'SRP231535', 
                        'SRP231534', 'SRP231533', 'SRP231129', 'SRP231128', 
                        'SRP231127', 'SRP231126', 'SRP231120', 'SRP231119',
                        "SRP231118", "SRP231117", 'SRP231035', 'SRP231034', 
                        'SRP231033', 'SRP231032', 'SRP231031', 'SRP231030', 
                        'SRP231029', 'SRP231028', 'SRP231027', 'SRP231026',
                        'SRP231025', 'SRP231024', 'SRP231023', 'SRP231022', 
                        'SRP231021', 'SRP231020', 'SRP231019', 'SRP231018', 
                        'SRP231017', 'SRP231016', 'SRP231015', 'SRP231014',
                        'SRP231013', 'SRP231012', 'SRP231011', 'SRP231010', 
                        'SRP231009', "SRP230142", "SRP230135", "SRP230133")


# read in sra table and parse  --------------------------------------------

sra <- read_csv("inputs/make_metadata/20220407_sra_results.csv", show_col_types = FALSE) %>%
  janitor::clean_names() %>%
  filter(library_source == "TRANSCRIPTOMIC") %>%
  filter(library_strategy %in% c("RNA-Seq", "OTHER")) %>%
  filter(library_selection %in% c("cDNA", "other", "RANDOM", "unspecified", "RT-PCR", "RANDOM PCR", "PCR")) %>%
  filter(study_accessions %in% c("SRP345944", "SRP305697")) %>% # remove broken study accessions (non-jgi)
  filter(study_accessions %in% jgi_sras_that_fail) # remove JGI submitted broken study accessions

study_accessions <- unique(sra$study_accession)
length(study_accessions)
head(study_accessions)
write.table(study_accessions, "inputs/make_metadata/20220407_srp_ids.txt", col.names = F, row.names = F, quote = F)
