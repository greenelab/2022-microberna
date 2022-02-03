# Compendia of bacterial and archaeal RNA sequencing data

## Getting started with this repository

This repository uses conda to manage software installations. 
You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/en/latest/miniconda.html).

```
conda env create --name microberna --file environment.yml
conda activate microberna
```

Once the base environment is installed, the analysis is automated with snakemake. 

```
snakemake -s best_genome_reference.snakefile -j 1 --use-conda --rerun-incomplete
```

If you're using an HPC such as slurm, snakemake can parallelize job submission and modulate resource usage (RAM, CPUs).

```
snakemake -s best_genome_reference.snakefile -j 16 --use-conda --rerun-incomplete --latency-wait 15 -
-resources mem_mb=200000 --cluster "sbatch -t 720 -J comp -p bmm -n 1 -N 1 -c {threads} --mem={resources.mem_mb}" -k
```

## Background and goals

The goal of this repository is to generate an automated, lightweight, and generalized pipeline to create a compendium of isolate bacterial and archaeal RNA-seq data.
In model organisms like human and mouse, compendia like [recount2](https://www.nature.com/articles/nbt.3838) (and [recount3](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02533-6)), [GTEx](https://gtexportal.org/home/), and [The Cancer Genome Atlas (TCGA)](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga) have been rich community resources that have allowed researchers to interface with and draw insights from gene expression data. 
Standard pre-processing enables integrated analysis across data sets derived from different studies, tissues, or organisms.
This repository aims to produce a compendium that integrates the thousands of publicly available isolate bacterial and archaeal RNA-seq data via uniform pre-processing. 

## Products

While the exact output formats will become more clear as the pre-processing pipeline is finalized, the pipeline aims to produce annotation-dependent and annotation-independent products.

**Annotation-dependent products**

+ Strain-level taxonomy profiles: we will use `sourmash gather` to determine the taxonomic profile of each isolate RNA-seq data set using genomes in the GTDB database (rs202) and the human transcriptome (GCF_000001405.39_GRCh38.p13) as references.
+ Uncounted (unmapped) reads: we will output a fastq file of reads that did not map to their reference.
+ Ortholog counts: reference transcriptomes, and therefore resultant counts, will be annotated using eggnog. We will using `tximport` to produce ortholog-summarized counts for each RNA-seq sample. Eggnog provides annotations for orthologs in COG, GO, KEGG, EC, CAZy, BiGG Reaction, and PFAM databases.

**Annotation-independent products**

+ FracMinHash sketches (*k* = 21, 31, scaled = 1000, abundance tracking): MinHash sketching systemically subsamples sequencing reads to allow rapid comparisons against other sketches (other transcriptomes, metatranscriptomes, or reference genomes). Abundance tracking captures differences in expression, but can also be used to k-mer trim rare (abundance = 1) k-mers that likely represent sequencing error and that may defalte similarity or containment estimates.

## Summary of RNA-seq samples

As of January 2022, there were 54,445<sup>1</sup> publicly available bacterial and archaeal RNA seq samples from 5,751<sup>2</sup> experiments on the Sequence Read Archive (SRA).
These samples were attributable to 1,722 distinct user-supplied organism names.
The median number of times an organism name was observed was 7.

Using text-based matching, we converted organism names in SRA metadata to species names in the GTDB taxonomy.
Through text-based matching alone, 3,274 RNA-seq samples derived from 344 organism names did not have a reference genome in the GTDB database.
2,524 GTDB species were represented in the 51,171 samples with matches, with a median of 29 observations per species.
The majority of samples were Bacterial (50,182 Bacterial, 989 Archaeal).

## Generating reference transcriptomes

A major challenge associated with building this compendium is the selection of reference sequences and the integration of annotations between reference sequences.
Here, we focus on generating a reference transcriptome, as we have elected to perform transcript quantification using Salmon because it is lightweight and scales to thousands of samples.

The reference transcriptome is an important component of transcriptome quantification.
If the reference does not contain all of the sequences that are contained within the sample, quantification will be incomplete (see [here](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0180904) and [here](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0145861)).
This is nicely illustrated in this figure by [Lisa K. Johnson](https://johnsolk.github.io/blog/) (original notebook [here](https://github.com/WhiteheadLab/RNAseq_17killifish/blob/master/notebooks/alignment_rates_plot.ipynb)). 
Johnson demonstrates that has Jaccard similarity decreases between an RNA-seq sample and the reference used to quantify it, mapping rates decrease. 
As such, it is important to have both a complete (all genes present) and similar (high percent identity) reference transcriptome with regards to each sample.

![](https://i.imgur.com/nSFIKHO.png)

Given the importance of having a representative and accurate transcriptome, one goal of this repository is to determine the best approach for reference transcriptome generation for any publicly available bacterial or archaeal RNA-seq sample.
Thus far, we have identified two approaches for reference transcriptome generation that we plan to test: a *single pangenome* approach and a *best genome* approach. 

The *single pangenome* approach will generate a single pangenome using all "reference" genomes for a species.
We have selected GTDB rs202 as our reference genome database, as it encompasses some unculturable organisms (unlike RefSeq), but is still quality controlled more extensively than e.g. GenBank. 
Annotated features within this pangenome (CDS, tRNA, rRNA, tmRNA, ncRNA) will be used to create a reference transcriptome that will be used to quantify all RNA-seq libraries from this species.
We suspect that this approach might be successful at capturing microbes that have mixes of genes (e.g. accessory elements) which have not previously be observed in conjunction in a single genome before; a [similar method](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0145861) was successful in the analysis of RNA-seq samples from *Staphylococcus aureus* strains without reference genomes 

The *best genome* approach will select the best genome and use it as the reference transcriptome for each RNA-seq library.
To achieve this, we will use sourmash `gather` to identify the reference genome in GTDB that best covers the k-mers in a given RNA-seq library.
Using this information, we will download this genome and identify, index, and annotate its open reading frames.
To harmonize annotations across data sets, we will rely on ortholog annotations produced by e.g. eggnog.
One benefit of this approach is the direct use of model organism transcriptomes (e.g. *Escherichia coli* k12, *Pseudomonas aeruginosa* PAO1 or PA14), which could be exported as sub-compendia and would better integrate with existing and well-known community annotation resources.

Prior to evaluation, we can foresee different challenges associated with each approach.

Challenges with *single pangenome* approach:

+ Inappropriately collapsing gene sequences with high percent identity and different functions
+ Decreased mapping rates for less-closely related strains (unclear if this would be an issue in the pangenome context)
    + This could potentially be overcome with amino acid mapping to pangenome (e.g. with paladin), but this would require intermediate BAM files, which are cumbersome.
+ Integration of model organism annotations.

Challenges with *best genome* approach:

+ Harmonization of genes/orthologs after quantification
    + This could be overcome with gene maps and tximport, especially given eggnog and its associated annotation databases.
+ Mixes of accessory genes within an RNA-seq library that have not been observed in conjunction in reference genomes before.
+ Decreased mapping rates for strains that are less-closely related to the best reference, e.g. decreased ANI independent of the observed accessory genes.

## Evaluation

### Comparison of reference transcriptome generation methods

### Toy compendia


#### *Faecalibacterium prausnitzii* A2-165

Three studies representing 29 samples have been submitted to the SRA for *F. prausnitzii* A2-165. 
While some of the samples are contaminated with human or other organisms, only strain A2-165 was detected among *F. prausnitzii* reads.
These RNA-seq samples are a good test case to determine whether the pangenome approach works as well as the best genome approach.
We selected three samples, each from a different study, and used these to investigate this question.

| SRA identifier | Annotated strain |
|----------------|------------------|
| ERX4307280     | A2-165  |
| SRX10245671    | A2-165  |
| SRX3847835     | A2-165  |

See results [here](notebooks/20220201_benchmark_salmon_vs_star_one_reference.ipynb)

#### *Bacillus* 

Five studies representing 21 samples have been submitted to the SRA and labelled as *Bacillus pumilus*. 
We selected one sample from each of these studies to test. 
Four of the samples were annotated by the submitters as *B. pumilus*, illustrating the necessity of determining the reference transcriptome by sequence comparison against a database. 
The difference in both the total number of *Bacillus* strains detected, as well as the identity of those strains (not shown), suggests that each transcriptome captures a slightly different strain, and that none of those strains are not represented in a single reference genome.
The mix of strains present makes these samples a good test case to determine whether the pangenome approach allows more reads to be recovered over the best genome approach.

| SRA identifier | Top annotated strain | Total *Bacillus* strains annotated |
|----------------|------------------|------------------------|
| SRX4378884     | Bacillus altitudinis strain=BA06 | 5 |
| SRX5179263     | Bacillus altitudinis strain=BA06 | 6 | 
| SRX5678458     | Bacillus altitudinis strain=BA06 | 8 |
| SRX7088415     | Bacillus altitudinis strain=BA06 | 5 | 

#### *Pseudomonas aerruginoas*

Recently, compendia of *P. aeruginosa* were developed by mapping RNA-seq samples against PAO1 or PA14 references using salmon. 
We identified 6 samples analyzed by this study to use as test cases for this compendium.

| SRA identifier | Annotated strain |
|----------------|------------------|
| ERX3558803     | PAO1 |
| SRX540112      | PAO1 |
| SRX589549      | PA14 |
| SRX4624095     | PA14 |
| SRX396878      | clinical isolate |
| SRX5123759     | clinical isolate |

*P. aeruginosa* is a good test case for this pipeline because:
1. It has many reference genomes (5,211 in GTDB), so we can test:
    a. How accurately we predict the correct reference genome
    b. containment of the RNA-seq sample in a single or many reference genomes 
2. PAO1 and PA14 are important community reference strains, so they'll be good test cases for integrating annotations at various levels.
    a. Do we collapse any important annotations by using ortholog-level annotations?
    b. Can we simultaneously take advantage of existing annotations and ortholog annotations?
3. The compendia referenced above have some nice controls, and we may be able to take advantage of some of those controls instead of recreating our own (e.g. count estimation with mapping vs. quasi-mapping).

### Impact of Salmon k-mer size on mapping rates and accuracy

Salmon in part uses k-mers in RNA-seq reads and transcripts in the reference transcriptome to quantify transcript abundance in a sample. 
One user-controllable parameter is the minimum k-size for an acceptable match.
We need to test how this parameter behaves in relation to the read size and the similarity between the RNA-seq sample and the selected reference transcriptome.

The salmon documentation [states](https://salmon.readthedocs.io/en/latest/salmon.html?highlight=k-mer#preparing-transcriptome-indices-mapping-based-mode):

> This will build the mapping-based index, using an auxiliary k-mer hash over k-mers of length 31. 
> While the mapping algorithms will make used of arbitrarily long matches between the query and reference, the k size selected here will act as the minimum acceptable length for a valid match. 
> Thus, a smaller value of k may slightly improve sensitivty. 
> We find that a k of 31 seems to work well for reads of 75bp or longer, but you might consider a smaller k if you plan to deal with shorter reads. 
> Also, a shoter value of k may improve sensitivity even more when using selective alignment (enabled via the –validateMappings flag). 
> So, if you are seeing a smaller mapping rate than you might expect, consider building the index with a slightly smaller k.

**Read length**: Estimating average read length by dividing SRA metadata Total Bases by Total Spots,

+ Only 387 samples have an average read length < 31, so a k size of 31 would be possible for the majority of data sets.
+ In contrast, 13,849 have an average read length < 75

It may be prudent to require read length > 31 base pairs after trimming for inclusion in the compendia. In either case, we need to test k size across read lengths to see if there are biases.

**ANI between genome in the RNA-seq sample and the reference transcriptome**: Decreased relatedness between the organism represented by the RNA-seq sample and that represented by the reference transcriptome should decrease mapping rates, but decreasing the minimum acceptable k-mer length could improve mapping rates in this case. We need to test whether decreasing the k size increases mapping rate without increasing off target mapping, for both RNA-seq samples that are closely and distantly related to the selected reference.

It would be great if there is a universal best k-size to use.

I will probably evaluate this after I have gather results for a few species. An ideal data set would have varied jaccard similarity relative to the reference pangenome (or other reference transcriptome used). 

## Potentially interesting use cases for such a compendia

**General ideas**

+ identification of conserved operons across species
+ conserved and divergent ortholog co-expression patterns across species
+ pathway enrichment analysis across environments (if most samples have metadata)
+ distribution of prophage, antibiotic resistance, gene expression
+ are there any contaminants, or contaminant-host pairs, that occur more frequently than would be expected by chance?

**Single pangenome ideas**

+ How many genes in the pangenome are observed in the transcriptome? How does this vary by species?
+ Estimate size, etc. of pangenomes for all species that are in GTDB and have gene expression data sets.
+ Identify horizontally transferred genes that are shared (and expressed) between species

## Footnotes
<sup>1</sup>Search criteria on the SRA ("Bacteria"[Organism] OR "Archaea"[Organism]) AND ("biomol rna"[Properties] AND "platform illumina"[Properties] AND "filetype fastq"[Properties]), and additionally filtered on Library Source "Transcriptomic", Library Strategy "RNA-Seq" or "OTHER", and Library Selection "cDNA", "other", "RANDOM", "unspecified", "RT-PCR", "RANDOM PCR", "PCR"  
<sup>2</sup>Measured as distinct study accessions
