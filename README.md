# Compendia of bacterial and archaeal RNA sequencing data

## Getting started with this repository

This repository uses conda to manage software installations. 
You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/en/latest/miniconda.html).

```
conda env create --name microberna --file environment.yml
conda activate refinebio
```

## Generating reference transcriptomes

The reference transcriptome is a quintessential piece for transcriptome quantification.
Tools like Salmon use exact k-mer matching between reads in a sample and transcripts in a reference transcriptome to quantify transcript abundances.
(Add notes from Lisa's killifish project about decreasing mapping rate with increasing phylogenetic distance from reference.)

Given the importance of having a representative and accurate transcriptome, one goal of this repository is to determine the best approach for reference transcriptome generation for any publicly available bacterial or archaeal RNA-seq data set.
Thus far, we have identified two approaches for reference transcriptome generation that we plan to test: a *single pangenome* approach and a *best genome* approach. 

The *single pangenome* approach will generate a single pangenome using all "reference" genomes for a species.
We have selected GTDB rs202 as our reference genome database, as it encompasses some unculturable organisms (unlike RefSeq), but is still quality controlled more extensively than e.g. GenBank. 
This genes within this pangenome will be used to create a reference transcriptome that will be used to quantify all RNA-seq libraries from this species.

The *best genome* approach will select the best genome and use it as the reference transcriptome for each RNA-seq library.
To achieve this, we will use sourmash `gather` to identify the reference genome in GTDB that best covers the k-mers in a given RNA-seq library.
Using this information, we will download this genome and identify, index, and annotate its open reading frames.
To harmonize annotations across data sets, we will rely on ortholog annotations in e.g. eggnog.

Prior to evaluation, we can foresee different challenges and benefits associated with each approach.

Challenges with *single pangenome* approach
+ inappropriately collapsing gene sequences with high percent identity and different functions
+ decreased mapping rates for less-closely related strains (unclear if this would be an issue in the pangenome context)
    + could be overcome with amino acid mapping to pangenome (e.g. with paladin), but this would require intermediate BAM files, which are cumbersome.
+ determining the best time to make gene:ortholog map for each genome. Should this just be done across all genomes in GTDB for a given species? Or should we wait and do this at the very end, when we know which genomes will actually be used as references? I'm inclined toward doing it at the beginning; once a species is observed in RNA-seq data, trigger the download of all GTDB genomes for that species and make a pangenome. Orrrr a pangenome may not be necessary. If we're using eggnog to annotate the orthologs, that should provide standardized names across genomes. 
+ But single pangenome approach could better capture microbes that have mixes of genes (e.g. accessory elements) which have not previously be observed in conjunction in a single genome before.

Challenges with *best genome* approach
+ harmonization of genes/orthologs after quantification
    + could be overcome with gene maps and tximport, especially given eggnog and its associated annotation databases.
+ Mixes of accessory genes within an RNA-seq library that have not been observed in conjunction in reference genomes before.

## Evaluation


## Potentially interesting use cases for such a compendia

**General ideas**

+ identification of conserved operons across species
+ pathway enrichment analysis across environments (if most samples have metadata)
+ distribution of prophage, antibiotic resistance, gene expression

**Single pangenome ideas**
+ how  many genes in the pangenome are observed in the transcriptome? How does this vary by species?
+ estimate size, etc. of pangenomes for all species that are in GTDB and have gene expression data sets.
