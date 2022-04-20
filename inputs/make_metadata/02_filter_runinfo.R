library(dplyr)
library(readr)

# this script reads in the ENA runinfo table, filters it back down to the
# transcriptome samples we are investigating, and writes it back out.
# It also investigates the size of the metadata table, and the identifiers
# that should be used for this project.


# read in and filter data -------------------------------------------------

# these samples will need to be filtered out, as they were not included in the 
# runinfo tables.
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

# read in sra table, filtering to transcriptome samples and removing samples
# that were not in the runinfo table
sra <- read_csv("inputs/make_metadata/20220407_sra_results.csv", show_col_types = FALSE) %>%
  janitor::clean_names() %>%
  filter(library_source == "TRANSCRIPTOMIC") %>%
  filter(library_strategy %in% c("RNA-Seq", "OTHER")) %>%
  filter(library_selection %in% c("cDNA", "other", "RANDOM", "unspecified", "RT-PCR", "RANDOM PCR", "PCR")) %>%
  filter(study_accessions %in% c("SRP345944", "SRP305697")) %>% # remove broken study accessions (non-jgi)
  filter(study_accessions %in% jgi_sras_that_fail) # remove JGI submitted broken study accessions

# read in run info and filter to ids in the sra table
runinfo <- read_tsv("inputs/make_metadata/20220407_srp_ids_all.runinfo.tsv")

runinfo_filt <- runinfo %>%
  filter(experiment_accession %in% sra$experiment_accession | secondary_sample_accession %in% sra$sample_accession) %>%
  distinct() %>%
  filter(library_source == "TRANSCRIPTOMIC") %>%
  filter(library_strategy %in% c("RNA-Seq", "OTHER")) %>%
  filter(library_selection %in% c("cDNA", "other", "RANDOM", "unspecified", "RT-PCR", "RANDOM PCR", "PCR")) %>%
  filter(!is.na(fastq_ftp))  %>% # some ftp links are missing, remove samples (removes 184 samples)
  filter(read_count >= 1000000)  %>% # filter to samples that have at least 1 million reads  (removes 2083 samples)
  mutate(read_length = base_count/read_count) %>% 
  filter(read_length > 31)  # filter very short reads (removes 350 samples)

# write out filtered metadata file
write_tsv(runinfo_filt, "inputs/20220407_runinfo.tsv") # 59239 samples

# look at filtered runinfo ------------------------------------------------

# look at composition, determine which identifier to use
dim(runinfo_filt)
dim(runinfo)
dim(sra)
length(unique(runinfo_filt$experiment_accession))
length(unique(runinfo_filt$run_accession))
length(unique(runinfo_filt$secondary_sample_accession))
length(unique(runinfo_filt$sample_accession))
length(unique(runinfo_filt$fastq_ftp))

tmp <- runinfo_filt %>%
  select(run_accession, fastq_ftp) %>%
  group_by(fastq_ftp) %>%
  tally() %>%
  arrange(desc(n))

# Secondary sample accession: some secondary sample accessions are repeated,
#                             but represent distinct sequencing libraires,
#                             as is clear from the run_alias (ex. SRS4756446)
# Experiment accession: some experiment accessions are repeated, but represent
#                       distinct sequencing libraries, as is clear by
#                       the metadata field run_alias (example SRX1090218).
# Run accessions: not repeated...currently best option

# which fastq_ftp links are missing?
tmp <- runinfo_filt %>%
  filter(is.na(fastq_ftp))


tmp <- runinfo_filt %>%
  select(experiment_accession, experiment_title, sample_title, read_count) %>%
  arrange(desc(read_count))
View(tmp)

library(ggplot2)
library(scales)
ggplot(runinfo_filt, aes(x = read_length)) +
  geom_histogram(bins = 50) +
  theme_minimal() +
  xlim(0, 100)

table(runinfo_filt$read_count > 500000)
