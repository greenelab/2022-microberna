import csv
import sys
import urllib.request
import pandas as pd

m = pd.read_csv("inputs/toy_metadata/toy_bp.tsv", sep = "\t", header = 0)
SRA = list(m['experiment_accession'])

TMPDIR = "/scratch/tereiter/" # TODO: update tmpdir based on computing env, or remove tmpdir invocation in resources
#GTDB_SPECIES = ['s__Faecalibacterium_prausnitzii_C']

gtdb_lineages = pd.read_csv("", header = 0)
GTDB_ACCS = 

class Checkpoint_RnaseqToReference:
    """
    Define a class a la genome-grist to simplify file specification
    from checkpoint (e.g. solve for {acc} wildcard). This approach
    is documented at this url:
    http://ivory.idyll.org/blog/2021-snakemake-checkpoints.html
    """
    def __init__(self, pattern, sra_samples):
        self.pattern = pattern
        self.sra_samples = sra_samples


    def get_sra_to_species(self, sra):
        sra_to_ref_csv = f'outputs/rnaseq_sourmash_gather_to_ref_species/{sra}.csv'
        assert os.path.exists(sra_to_ref_csv)

        # there should only be one sra_to_ref, as this reads in the csv
        # that records the best species-level match for each SRA sample
        with open(sra_to_ref_csv, 'rt') as fp:
           r = csv.DictReader(fp)
           for row in r:
               sra_to_ref = row['species_no_space']
               break            # maybe assert there is only one line?

        return sra_to_ref


    def __call__(self, w):
        global checkpoints

        # wait for the results of 'query_to_species_db';
        # this will trigger exception until that rule has been run.
        # look at do_sample in genome_grist Snakefile for cleaner

        # this will wait for all SRA samples to be done... which is possibly
        # unnecessary...
        for sra in self.sra_samples:
            d = dict(sra=sra)
            ww = snakemake.io.Wildcards(fromdict=d)

            checkpoints.rnaseq_sample_select_best_species_reference.get(**ww)

        # parse accessions in gather output file
        #ref_to_sra_res = self.get_sra_to_species(w.sra)

        # what you want here is essentially a double for loop -
        # for each sra,
        #    get the species that matches to that SRA
        #        and then shove that into the pattern
        results = []
        for sra in self.sra_samples:
            found = False
            for p in self.pattern:
                if sra in p:
                    assert not found
                    # this is the one and only element of the pattern list
                    # where we are going to put this 'gtdb_species' into the
                    # pattern.
                    gtdb_species = self.get_sra_to_species(sra)
                    p = expand(p, gtdb_species=gtdb_species)
                    results.append(p)
                    found = True
        #p = expand(self.pattern, gtdb_species=ref_to_sra_res, **w)

        #print('ZZZ', results)
        
        results = sum(results, []) # turn list of lists into a list of character strings
        return results


class Checkpoint_GrabAccessions:
    """
    Define a class to simplify file specification from checkpoint 
    (e.g. solve for {acc} wildcard without needing to specify a function
    for each arm of the DAG that uses the acc wildcard). 
    This approach is documented at this url: 
    http://ivory.idyll.org/blog/2021-snakemake-checkpoints.html
    """
    def __init__(self, pattern):
        self.pattern = pattern

    def get_genome_accs(self, gtdb_species):
        acc_csv = f'outputs/gtdb_genomes_by_species/{gtdb_species}.csv'
        assert os.path.exists(acc_csv)

        genome_accs = []
        with open(acc_csv, 'rt') as fp:
           r = csv.DictReader(fp)
           for row in r:
               acc = row['ident']
               genome_accs.append(acc)
        print(f'loaded {len(genome_accs)} accessions from {acc_csv}.')

        return genome_accs

    def __call__(self, w):
        global checkpoints

        # wait for the results of rule 'grab_species_accessions'; 
        # this will trigger exception until that rule has been run.
        checkpoints.grab_species_accessions.get(**w)

        # parse accessions in gather output file
        genome_accs = self.get_genome_accs(w.gtdb_species)

        p = expand(self.pattern, acc=genome_accs, **w)
        return p


rule all:
    input:
        # generate pangenome and annotate with eggnog
        #expand("outputs/gtdb_genomes_roary_eggnog/{gtdb_species}.emapper.annotations", gtdb_species = GTDB_SPECIES),
        # generate pantranscriptome and index it with salmon
        #expand("outputs/gtdb_genomes_salmon_index/{gtdb_species}/info.json", gtdb_species = GTDB_SPECIES)
        # gather RNAseq sample
        #expand("outputs/rnaseq_sourmash_gather/{sra}_gtdb_k31.csv", sra = SRA),
        Checkpoint_RnaseqToReference(expand("outputs/rnaseq_salmon/{{gtdb_species}}/{sra}_quant/quant.sf", sra = SRA), SRA),

##############################################################
## Generate reference transcriptome using pangenome analysis
##############################################################

# TODO: update to GTDB RS207 taxonomy sheet
# probably need to package with git repo as an input to expand over uncompressed directory
rule download_gtdb_lineages:
    output: "inputs/gtdb-rs202.taxonomy.v2.csv"
    resources:
        mem_mb = 4000
    threads: 1
    shell:"""
    wget -O {output} https://osf.io/p6z3w/download
    """
 
# The outputs of this rule are used by the class Checkpoint_GrabAccessions.
# The new wildcard, acc, is encoded in the "ident" column of the output
# file "accessions". 
checkpoint grab_species_accessions:
    input: gtdb_lineages = "inputs/gtdb-rs202.taxonomy.v2.csv"
    output: accessions = "outputs/gtdb_genomes_by_species/{gtdb_species}.csv"
    resources:
        mem_mb = 4000
    threads: 1
    conda: "envs/tidyverse.yml"
    script: "scripts/grab_species_accessions.R"


rule bakta_download_db:
    output: "inputs/bakta_db/db/version.json"
    threads: 1
    resources: mem_mb = 4000
    params: outdir = "inputs/bakta_db"
    conda: "envs/bakta.yml"
    shell:'''
    bakta_db download --output {params.outdir}
    '''

rule download_charcoal_gtdb_rs207_genomes:
    output: "inputs/charcoal_gtdb_rs207_genomes.tar.gz"
    threads: 1
    resources: mem_mb = 1000
    shell:'''
    wget -O {output} URL # TODO: ADD URL
    '''

rule decompress_charcoal_gtdb_rs207_genomes:
    input: "inputs/charcoal_gtdb_rs207_genomes.tar.gz"
    output: expand("inputs/charcoal_gtdb_rs207_genomes/{all_accs}_genomic.fna.gz.clean.fa.gz", all_accs = ALL_ACCS)
    params: outdir = "outputs/charcoal_gtdb_rs207_genomes/"
    resources: mem_mb = 1000
    shell:'''
    tar xf {input} -C {params.outdir}
    '''

rule bakta_annotate_gtdb_genomes:
    input: 
        fna=ancient("inputs/charcoal_gtdb_rs207_genomes/{acc}_genomic.fna.gz.clean.fa.gz"),
        db="inputs/bakta_db/db/version.json",
    output: 
        "outputs/gtdb_genomes_bakta/{gtdb_species}/{acc}.faa",
        "outputs/gtdb_genomes_bakta/{gtdb_species}/{acc}.gff3",
        "outputs/gtdb_genomes_bakta/{gtdb_species}/{acc}.fna",
        "outputs/gtdb_genomes_bakta/{gtdb_species}/{acc}.ffn",
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 16000 ,
        tmpdir= TMPDIR
    benchmark: "benchmarks/bakta_{gtdb_species}/{acc}.txt"
    conda: 'envs/bakta.yml'
    params: 
        dbdir="inputs/bakta_db/db/",
        outdir = lambda wildcards: 'outputs/gtdb_genomes_bakta/' + wildcards.gtdb_species,
    threads: 1
    shell:'''
    bakta --threads {threads} --db {params.dbdir} --prefix {wildcards.acc} --output {params.outdir} \
        --locus {wildcards.acc} --locus-tag {wildcards.acc} --keep-contig-headers {input.fna}
    '''

rule cat_annotated_sequences:
    input: ancient(Checkpoint_GrabAccessions("outputs/gtdb_genomes_bakta/{{gtdb_species}}/{acc}.ffn"))
    output: "outputs/gtdb_genomes_annotated_comb/{gtdb_species}_all_annotated_seqs.fa"
    threads: 1
    resources: 
        mem_mb=2000,
        tmpdir=TMPDIR
    benchmark:"benchmarks/cat_annotated_seqs_{gtdb_species}.txt"
    shell:'''
    cat {input} > {output}
    '''

# note clustering might not be necessary because of the pufferfish index underlying salmon indexing. 
# but in any case it should help with downstream accounting, so keep.
rule cluster_annotated_sequences:
    input: "outputs/gtdb_genomes_annotated_comb/{gtdb_species}_all_annotated_seqs.fa"
    output: "outputs/gtdb_genomes_annotated_comb/{gtdb_species}_clustered_annotated_seqs.fa"
    threads: 1
    resources: 
        mem_mb=16000,
        tmpdir=TMPDIR
    conda: "envs/cdhit.yml"
    benchmark:"benchmarks/cluster_annotated_{gtdb_species}.txt"
    shell:'''
    cd-hit-est -c .95 -d 0 -i {input} -o {output}
    '''

rule translate_clustered_sequences_for_annotations:
    input: "outputs/gtdb_genomes_annotated_comb/{gtdb_species}_clustered_annotated_seqs.fa" 
    output: 'outputs/gtdb_genomes_annotated_comb/{gtdb_species}_clustered_annotated_seqs.faa'
    conda: 'envs/emboss.yml'
    resources:
        mem_mb = 16000,
        tmpdir = TMPDIR
    benchmark: "benchmarks/transeq_{gtdb_species}.txt"
    threads: 2
    shell:'''
    transeq {input} {output}
    '''

rule eggnog_download_db:
    output: "inputs/eggnog_db/eggnog.db"
    threads: 1   
    resources: 
        mem_mb = 4000,
        tmpdir=TMPDIR
    params: datadir = "inputs/eggnog_db"
    conda: "envs/eggnog.yml"
    shell:'''
    download_eggnog_data.py -H -d 2 -y --data_dir {params.datadir}
    '''

rule eggnog_annotate_clustered_sequences:
    input: 
        faa = "outputs/gtdb_genomes_annotated_comb/{gtdb_species}_clustered_annotated_seqs.faa",
        db = 'inputs/eggnog_db/eggnog.db'
    output: "outputs/gtdb_genomes_annotated_comb_eggnog/{gtdb_species}/{gtdb_species}.emapper.annotations"
    conda: 'envs/eggnog.yml'
    resources:
        mem_mb = 32000,
        tmpdir=TMPDIR
    benchmark: "benchmarks/eggnog_{gtdb_species}.txt"
    threads: 8
    params: 
        outdir = lambda wildcards: "outputs/gtdb_genomes_annotated_comb_eggnog/" + wildcards.gtdb_species,
        dbdir = "inputs/eggnog_db"
    shell:'''
    mkdir -p tmp/
    emapper.py --cpu {threads} -i {input.faa} --output {wildcards.gtdb_species} \
       --output_dir {params.outdir} -m hmmer -d none --tax_scope auto \
       --go_evidence non-electronic --target_orthologs all --seed_ortholog_evalue 0.001 \
       --seed_ortholog_score 60 --override --temp_dir tmp/ \
       -d 2 --data_dir {params.dbdir}
    '''

#########################################################
## generate pangenome
#########################################################

rule roary_determine_pangenome:
    input: ancient(Checkpoint_GrabAccessions("outputs/gtdb_genomes_bakta/{{gtdb_species}}/{acc}.gff3"))
    output: "outputs/gtdb_genomes_roary/{gtdb_species}/pan_genome_reference.fa"
    conda: 'envs/roary.yml'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000,
        tmpdir=TMPDIR
    benchmark: "benchmarks/roary_{gtdb_species}.txt"
    threads: 2
    params: 
        outdir = lambda wildcards: 'outputs/gtdb_genomes_roary/' + wildcards.gtdb_species
    shell:'''
    roary -r -e -n -f {params.outdir} -p {threads} -z {input}
    '''


############################################################
## Generate decoy sequences
############################################################

rule define_genome_regions_as_bed:
    input: ancient("outputs/gtdb_genomes_bakta/{gtdb_species}/{acc}.fna"),
    output: "outputs/gtdb_genomes_regions_bed/{gtdb_species}/{acc}_region.sizes"
    threads: 1
    resources: 
        mem_mb=2000,
        tmpdir=TMPDIR
    conda: "envs/samtools.yml"
    benchmark:"benchmarks/genome_regions_{gtdb_species}/{acc}.txt"
    shell:'''    
    samtools faidx {input}
    cut -f 1,2 {input}.fai | sort -n -k1 > {output}
    '''

rule filter_bakta_gff:
    input:  gff=ancient("outputs/gtdb_genomes_bakta/{gtdb_species}/{acc}.gff3"),
    output: gff=temp("outputs/gtdb_genomes_gff/{gtdb_species}/{acc}_filtered.gff3"),
    threads: 1
    resources: 
        mem_mb=2000,
        tmpdir=TMPDIR
    conda: "envs/rtracklayer.yml"
    benchmark:"benchmarks/filter_gff_{gtdb_species}/{acc}.txt"
    script: "scripts/filter_gff.R"

rule sort_bakta_gff:
    input: gff="outputs/gtdb_genomes_gff/{gtdb_species}/{acc}_filtered.gff3",
    output: gff="outputs/gtdb_genomes_gff/{gtdb_species}/{acc}_filtered_sorted.gff3",
    threads: 1
    resources: 
        mem_mb=2000,
        tmpdir=TMPDIR
    conda: "envs/bedtools.yml"
    benchmark:"benchmarks/sort_gff_{gtdb_species}/{acc}.txt"
    shell:'''
    sortBed -i {input} > {output}
    '''

rule complement_gff:
    input:
        gff="outputs/gtdb_genomes_gff/{gtdb_species}/{acc}_filtered_sorted.gff3",
        sizes= "outputs/gtdb_genomes_regions_bed/{gtdb_species}/{acc}_region.sizes"
    output:"outputs/gtdb_genomes_intergenic_bed/{gtdb_species}/{acc}.bed"
    threads: 1
    resources: 
        mem_mb=2000,
        tmpdir=TMPDIR
    conda: "envs/bedtools.yml"
    benchmark:"benchmarks/complement_gff_{gtdb_species}/{acc}.txt"
    shell:'''
    bedtools complement -i {input.gff} -g {input.sizes} > {output}
    '''

rule extract_intergenic_from_genome:
    input:
        fna=ancient("outputs/gtdb_genomes_bakta/{gtdb_species}/{acc}.fna"),
        bed="outputs/gtdb_genomes_intergenic_bed/{gtdb_species}/{acc}.bed"
    output:"outputs/gtdb_genomes_intergenic_seqs/{gtdb_species}/{acc}.fa"
    threads: 1
    resources: 
        mem_mb=2000,
        tmpdir=TMPDIR
    conda: "envs/bedtools.yml"
    benchmark:"benchmarks/extract_intergenic_{gtdb_species}/{acc}.txt"
    shell:'''
    bedtools getfasta -fi {input.fna} -bed {input.bed} -name -s > {output}
    '''

rule cat_intergenic_sequences:
    input: ancient(Checkpoint_GrabAccessions("outputs/gtdb_genomes_intergenic_seqs/{{gtdb_species}}/{acc}.fa")),
    output: "outputs/gtdb_genomes_intergenic_comb/{gtdb_species}_all_intergenic_seqs.fa"
    threads: 1
    resources: 
        mem_mb=2000,
        tmpdir=TMPDIR
    benchmark:"benchmarks/cat_intergenic_{gtdb_species}.txt"
    shell:'''
    cat {input} > {output}
    '''

rule cluster_intergenic_sequences:
    input: "outputs/gtdb_genomes_intergenic_comb/{gtdb_species}_all_intergenic_seqs.fa"
    output: "outputs/gtdb_genomes_intergenic_comb/{gtdb_species}_clustered_intergenic_seqs.fa"
    threads: 1
    resources: 
        mem_mb=16000,
        tmpdir=TMPDIR
    conda: "envs/cdhit.yml"
    benchmark:"benchmarks/cluster_intergenic_{gtdb_species}.txt"
    shell:'''
    cd-hit-est -c 1 -i {input} -o {output}
    '''

###################################################
## Generate decoy-aware salmon transcriptome index
###################################################

rule grab_intergenic_sequence_names:
    input: "outputs/gtdb_genomes_intergenic_comb/{gtdb_species}_clustered_intergenic_seqs.fa"
    output: "outputs/gtdb_genomes_intergenic_comb/{gtdb_species}_clustered_intergenic_seq_names.txt"
    threads: 1
    resources: 
        mem_mb=2000,
        tmpdir=TMPDIR
    benchmark:"benchmarks/grab_intergenic_names_{gtdb_species}.txt"
    shell:"""
    grep '>' {input} | awk 'sub(/^>/, "")' > {output}
    """

rule combine_pangenome_and_intergenic_sequences:
    input: 
        "outputs/gtdb_genomes_annotated_comb/{gtdb_species}_clustered_annotated_seqs.fa",
        "outputs/gtdb_genomes_intergenic_comb/{gtdb_species}_clustered_intergenic_seqs.fa"
    output: "outputs/gtdb_genomes_salmon_ref/{gtdb_species}.fa"
    threads: 1
    resources: 
        mem_mb=2000,
        tmpdir=TMPDIR
    benchmark:"benchmarks/combine_sequences_{gtdb_species}.txt"
    shell:'''
    cat {input} > {output}
    '''

rule index_transcriptome:
    input: 
        seqs=ancient("outputs/gtdb_genomes_salmon_ref/{gtdb_species}.fa"),
        decoys=ancient("outputs/gtdb_genomes_intergenic_comb/{gtdb_species}_clustered_intergenic_seq_names.txt")
    output: "outputs/gtdb_genomes_salmon_index/{gtdb_species}/info.json"
    threads: 1
    params: index_dir = lambda wildcards: "outputs/gtdb_genomes_salmon_index/" + wildcards.gtdb_species
    resources: 
        mem_mb=16000,
        tmpdir=TMPDIR
    conda: "envs/salmon.yml"
    benchmark:"benchmarks/salmon_index_{gtdb_species}.txt"
    shell:'''
    salmon index -t {input.seqs} -i {params.index_dir} --decoys {input.decoys} -k 31
    '''

#################################################################
## Download, process, and taxonomically annotate RNAseq samples
#################################################################

rule rnaseq_sample_download:
    """
    Need to figure out better download specification.
    I went to the ENA, searched six samples I cared about, and downloaded the
    metadata files, which contain the URLs. I assume there is a programmatic way
    to generate ENA metadata tables for accessions of interest. Will investigate
    more later.
    """
    output:
        reads="outputs/rnaseq_fastp/{sra}.fq.gz",
        json = "outputs/rnaseq_fastp/{sra}.fastp.json",
        html = "outputs/rnaseq_fastp/{sra}.fastp.html"
    params: tmp_base = lambda wildcards: "inputs/tmp_raw/" + wildcards.sra
    threads: 1
    resources:
        mem_mb=8000
    run:
        row = m.loc[m['experiment_accession'] == wildcards.sra]
        fastqs = row['fastq_ftp'].values[0]
        fastqs = fastqs.split(";")
        if len(fastqs) == 1:
            # single end data; download and stream directly to fastp for trimming.
            fastq = fastqs[0]
            shell("mkdir -p inputs/tmp_raw")
            if not os.path.exists(params.tmp_base + ".fastq.gz"):
                shell("wget -O {params.tmp_base}.fastq.gz ftp://{fastq}")

            shell("fastp -i {params.tmp_base}.fastq.gz --json {output.json} --html {output.html} -R {wildcards.sra} --stdout | gzip > {output.reads}")

            # check that the file exists, and if it does, remove raw fastq files
            if os.path.exists(output.reads):
                os.remove(params.tmp_base + ".fastq.gz")

        else:
            # paired end data; download both files, interleave, and then remove files
            fastq_1 = fastqs[0]
            fastq_2 = fastqs[1]
            shell("mkdir -p inputs/tmp_raw")
            if not os.path.exists(params.tmp_base + "_1.fastq.gz"):
                shell("wget -O {params.tmp_base}_1.fastq.gz ftp://{fastq_1}")

            if not os.path.exists(params.tmp_base + "_2.fastq.gz"):
                shell("wget -O {params.tmp_base}_2.fastq.gz ftp://{fastq_2}")

            shell("fastp -i {params.tmp_base}_1.fastq.gz -I {params.tmp_base}_2.fastq.gz --json {output.json} --html {output.html} -R {wildcards.sra} --stdout | gzip > {output.reads}")

            # check that the file exists, and if it does, remove raw fastq files
            if os.path.exists(output.reads):
                os.remove(params.tmp_base + "_1.fastq.gz")
                os.remove(params.tmp_base + "_2.fastq.gz")
                

rule rnaseq_sample_sourmash_sketch:
    input: "outputs/rnaseq_fastp/{sra}.fq.gz",
    output: "outputs/rnaseq_sourmash_sketch/{sra}.sig"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000 ,
        tmpdir= TMPDIR
    threads: 1
    benchmark: "benchmarks/rnaseq/sourmash_sketch_{sra}.txt"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sketch dna -p k=21,k=31,k=51,scaled=1000,abund -o {output} --name {wildcards.sra} {input}
    '''

# update download link to RS207 when on OSF
rule sourmash_download_gtdb_database:
    output: "inputs/sourmash_dbs/gtdb-rs207.genomic.k31.zip" 
    resources:
        mem_mb = 1000,
        tmpdir = TMPDIR
    threads: 1
    shell:'''
    wget -O {output} https://osf.io/94mzh/download
    '''

rule sourmash_download_human_sig:
    output: "inputs/sourmash_dbs/GCF_000001405.39_GRCh38.p13_rna.sig"
    resources:
        mem_mb = 1000,
        tmpdir = TMPDIR
    threads: 1
    shell:'''
    wget -O {output} https://osf.io/anj6b/download
    '''

rule rnaseq_sample_sourmash_gather_against_gtdb:
    input:
        sig="outputs/rnaseq_sourmash_sketch/{sra}.sig",
        db="inputs/sourmash_dbs/gtdb-rs207.genomic.k31.zip",
        human="inputs/sourmash_dbs/GCF_000001405.39_GRCh38.p13_rna.sig"
    output: "outputs/rnaseq_sourmash_gather/{sra}_gtdb_k31.csv"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000 ,
        tmpdir= TMPDIR
    threads: 1
    benchmark: "benchmarks/rnaseq/sourmash_gather_k31_{sra}.txt"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash gather -o {output} --scaled 2000 -k 31 {input.sig} {input.db} {input.human}
    '''

checkpoint rnaseq_sample_select_best_species_reference:
    """
    Corresponds to class Checkpoint_RnaseqToReference at top of file.
    This will read in the gather results, and output one dataframe per RNA seq sample. 
    The dataframe will have the SRA accession and the reference genome accession.
    wildcard {sra} will already be defined as a list of RNAseq SRA accessions at the
    beginning of the Snakefile, and the checkpoint class will output an amalgamation 
    of wildcards {gtdb_species}-{sra}, which will be called {gtdb_species_to_sra} in 
    the rule all, but the hyphen should allow for the wildcards to be solved separately.
    Example:
    https://github.com/taylorreiter/2021-metapangenome-example/blob/main/Snakefile#L90
    """
    input:
        gtdb_lineages="inputs/gtdb-rs202.taxonomy.v2.csv",
        gather="outputs/rnaseq_sourmash_gather/{sra}_gtdb_k31.csv"
    output: sra_to_ref_species="outputs/rnaseq_sourmash_gather_to_ref_species/{sra}.csv"
    conda: "envs/tidyverse.yml"
    resources:
        mem_mb = 2000,
        tmpdir= TMPDIR
    threads: 1
    benchmark: "benchmarks/rnaseq/select_species_genome_{sra}.txt"
    script: "scripts/select_best_species_reference.R"

rule rnaseq_quantify_against_species_pangenome:
    input: 
        sra_to_ref_species="outputs/rnaseq_sourmash_gather_to_ref_species/{sra}.csv",
        index = ancient("outputs/gtdb_genomes_salmon_index/{gtdb_species}/info.json"),
        reads = "outputs/rnaseq_fastp/{sra}.fq.gz",
        eggnog = "outputs/gtdb_genomes_annotated_comb_eggnog/{gtdb_species}/{gtdb_species}.emapper.annotations" # not actually needed here, but trick snakemake into making these annots
    output: "outputs/rnaseq_salmon/{gtdb_species}/{sra}_quant/quant.sf"
    params: 
        index_dir = lambda wildcards: "outputs/gtdb_genomes_salmon_index/" + wildcards.gtdb_species,
        out_dir = lambda wildcards: "outputs/rnaseq_salmon/" + wildcards.gtdb_species + "/" + wildcards.sra + "_quant" 
    conda: "envs/salmon.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000 ,
        tmpdir= TMPDIR
    threads: 1
    benchmark: "benchmarks/rnaseq/salmon_quantify_{gtdb_species}-{sra}.txt"
    shell:'''
    salmon quant -i {params.index_dir} -l A -r {input.reads} -o {params.out_dir} --validateMappings --writeUnmappedNames
    '''
