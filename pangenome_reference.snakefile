import csv
import sys
import urllib.request
import pandas as pd

m = pd.read_csv("inputs/toy_metadata/toy_fpc.tsv", sep = "\t", header = 0)
SRA = m['experiment_accession']

TMPDIR = "/scratch/tereiter/" # TODO: update tmpdir based on computing env, or remove tmpdir invocation in resources
#GTDB_SPECIES = ['s__Pseudomonas_aeruginosa']
#GTDB_SPECIES = ['s__Bacillus_pumilus']
GTDB_SPECIES = ['s__Faecalibacterium_prausnitzii_C']
class Checkpoint_RnaseqToReference:
    """
    Define a class a la genome-grist to simplify file specification
    from checkpoint (e.g. solve for {acc} wildcard). This approach
    is documented at this url:
    http://ivory.idyll.org/blog/2021-snakemake-checkpoints.html
    """
    def __init__(self, pattern):
        self.pattern = pattern

    def get_acc_dbs(self, sra):
        sra_to_ref_csv = f'outputs/rnaseq_sourmash_gather_to_ref_species/{sra}.csv'
        assert os.path.exists(sra_to_ref_csv)

        # there should only be one sra_to_ref, as this reads in the csv
        # that records the best species-level match for each SRA sample
        with open(sra_to_ref_csv, 'rt') as fp:
           r = csv.DictReader(fp)
           for row in r:
               sra_to_ref = row['sra_to_ref_species']

        return sra_to_ref

    def __call__(self, w):
        global checkpoints

        # wait for the results of 'query_to_species_db';
        # this will trigger exception until that rule has been run.
        checkpoints.rnaseq_sample_select_best_species_reference.get(**w)

        # parse accessions in gather output file
        ref_to_sra_res = self.rnaseq_sample_select_best_species_reference(w.sra)

        p = expand(self.pattern, gtdb_species_to_sra=ref_to_sra_res, **w)
        return p


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
        expand("outputs/gtdb_genomes_salmon_index/{gtdb_species}/info.json", gtdb_species = GTDB_SPECIES)
        # gather RNAseq sample
        #expand("outputs/rnaseq_sourmash_gather/{sra}_gtdb_k31.csv", sra = SRA),
        #Checkpoint_RnaseqToReference(expand("outputs/rnaseq_salmon/{sra}/{{gtdb_species_to_sra}}_quant/quant.sf", sra = SRA))

##############################################################
## Generate reference transcriptome using pangenome analysis
##############################################################

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
    output: 
        accessions = "outputs/gtdb_genomes_by_species/{gtdb_species}.csv",
        charcoal_lineages = "outputs/gtdb_genomes_by_species/{gtdb_species}_charcoal_lineages.csv"
    resources:
        mem_mb = 4000
    threads: 1
    conda: "envs/tidyverse.yml"
    script: "scripts/grab_species_accessions.R"

rule make_genome_info_csv:
    output:
        csv = 'inputs/gtdb_genomes/{gtdb_species}/{acc}.info.csv'
    conda: "envs/download_genomes.yml"
    resources:
        mem_mb = 8000
    threads: 1
    shell: """
    python scripts/genbank_genomes.py {wildcards.acc} --output {output.csv}
    """

rule download_matching_genomes:
    input:
        csvfile = ancient('inputs/gtdb_genomes/{gtdb_species}/{acc}.info.csv')
    output:
        genome = "inputs/gtdb_genomes/{gtdb_species}/{acc}_genomic.fna.gz"
    resources:
        mem_mb = 500
    threads: 1
    run:
        with open(input.csvfile, 'rt') as infp:
            r = csv.DictReader(infp)
            rows = list(r)
            assert len(rows) == 1
            row = rows[0]
            acc = row['acc']
            assert wildcards.acc.startswith(acc)
            url = row['genome_url']
            name = row['ncbi_tax_name']

            print(f"downloading genome for acc {acc}/{name} from NCBI...",
                file=sys.stderr)
            with open(output.genome, 'wb') as outfp:
                with urllib.request.urlopen(url) as response:
                    content = response.read()
                    outfp.write(content)
                    print(f"...wrote {len(content)} bytes to {output.genome}",
                        file=sys.stderr)


rule sourmash_download_gtdb_database:
    output: "inputs/sourmash_dbs/gtdb-rs202.genomic.k31.zip" 
    resources:
        mem_mb = 1000,
        tmpdir = TMPDIR
    threads: 1
    shell:'''
    wget -O {output} https://osf.io/94mzh/download
    '''

rule charcoal_generate_genome_list:
    input:  ancient(Checkpoint_GrabAccessions("inputs/gtdb_genomes/{{gtdb_species}}/{acc}_genomic.fna.gz"))
    output: "outputs/charcoal_conf/charcoal_{gtdb_species}.genome-list.txt"
    threads: 1
    resources:
        mem_mb=500,
        tmpdir = TMPDIR
    params: genome_dir = lambda wildcards: "inputs/gtdb_genomes/" + wildcards.gtdb_species
    shell:'''
    ls {params.genome_dir}/*fna.gz | xargs -n 1 basename > {output} 
    '''

rule charcoal_generate_conf_file:
    input:
        genomes = ancient(Checkpoint_GrabAccessions("inputs/gtdb_genomes/{{gtdb_species}}/{acc}_genomic.fna.gz")),
        genome_list = "outputs/charcoal_conf/charcoal_{gtdb_species}.genome-list.txt",
        charcoal_lineages = "outputs/gtdb_genomes_by_species/{gtdb_species}_charcoal_lineages.csv",
        db="inputs/sourmash_dbs/gtdb-rs202.genomic.k31.zip",
        db_lineages="inputs/gtdb-rs202.taxonomy.v2.csv"
    output:
        conf = "conf/charcoal-conf-{gtdb_species}.yml",
    params: 
        genome_dir = lambda wildcards: "inputs/gtdb_genomes/" + wildcards.gtdb_species,
        output_dir = lambda wildcards: "outputs/gtdb_genomes_charcoal/" + wildcards.gtdb_species 
    resources:
        mem_mb = 500,
        tmpdir = TMPDIR
    threads: 1
    run:
        with open(output.conf, 'wt') as fp:
            print(f"""\
output_dir: {params.output_dir}
genome_list: {input.genome_list}
genome_dir: {params.genome_dir}
provided_lineages: {input.charcoal_lineages}
match_rank: order
gather_db:
 - {input.db} 
lineages_csv: {input.db_lineages} 
strict: 1
""", file=fp)


checkpoint charcoal_decontaminate_genomes:
    input:
        genomes = ancient(Checkpoint_GrabAccessions("inputs/gtdb_genomes/{{gtdb_species}}/{acc}_genomic.fna.gz")),
        genome_list = "outputs/charcoal_conf/charcoal_{gtdb_species}.genome-list.txt",
        conf = "conf/charcoal-conf-{gtdb_species}.yml",
        charcoal_lineages = "outputs/gtdb_genomes_by_species/{gtdb_species}_charcoal_lineages.csv",
        db="inputs/sourmash_dbs/gtdb-rs202.genomic.k31.zip",
        db_lineages="inputs/gtdb-rs202.taxonomy.v2.csv"
    output: directory("outputs/gtdb_genomes_charcoal/{gtdb_species}/"),
    resources:
        mem_mb = 275000,
        tmpdir = TMPDIR
    benchmark: "benchmarks/charcoal_{gtdb_species}.txt"
    threads: 32
    conda: "envs/charcoal.yml"
    shell:'''
    python -m charcoal run {input.conf} -j {threads} clean --nolock --latency-wait 15 --rerun-incomplete
    '''

rule bakta_download_db:
    output: "inputs/bakta_db/db/version.json"
    threads: 1
    resources: mem_mb = 4000
    params: outdir = "inputs/bakta_db"
    conda: "envs/bakta.yml"
    shell:'''
    bakta_db download --output {params.outdir}
    '''

rule bakta_annotate_gtdb_genomes:
    # TODO: change locus_tag to accept wildcards.acc if bakta/#92 gets implemented.
    input: 
        fna=ancient("outputs/gtdb_genomes_charcoal/{gtdb_species}/{acc}_genomic.fna.gz.clean.fa.gz"),
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

def checkpoint_charcoal_decontaminate_genomes1(wildcards):
    # checkpoint_output encodes the output dir from the checkpoint rule.
    checkpoint_output = checkpoints.charcoal_decontaminate_genomes.get(**wildcards).output[0]    
    file_names = expand("outputs/gtdb_genomes_bakta/{{gtdb_species}}/{acc}.ffn",
                        acc = glob_wildcards(os.path.join(checkpoint_output, "{acc}_genomic.fna.gz.clean.fa.gz")).acc)
    return file_names

rule cat_annotated_sequences:
    input: checkpoint_charcoal_decontaminate_genomes1
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
    cd-hit-est -c .95 -i {input} -o {output}
    '''

#########################################################
## generate pangenome
#########################################################

def checkpoint_charcoal_decontaminate_genomes2(wildcards):
    # checkpoint_output encodes the output dir from the checkpoint rule.
    checkpoint_output = checkpoints.charcoal_decontaminate_genomes.get(**wildcards).output[0]    
    file_names = expand("outputs/gtdb_genomes_bakta/{{gtdb_species}}/{acc}.gff3",
                        acc = glob_wildcards(os.path.join(checkpoint_output, "{acc}_genomic.fna.gz.clean.fa.gz")).acc)
    return file_names

rule roary_determine_pangenome:
    input: checkpoint_charcoal_decontaminate_genomes2,
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

rule translate_pangenome_for_annotations:
    input: 'outputs/gtdb_genomes_roary/{gtdb_species}/pan_genome_reference.fa' 
    output: 'outputs/gtdb_genomes_roary/{gtdb_species}/pan_genome_reference.faa' 
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

rule eggnog_annotate_pangenome:
    input: 
        faa = "outputs/gtdb_genomes_roary/{gtdb_species}/pan_genome_reference.faa",
        db = 'inputs/eggnog_db/eggnog.db'
    output: "outputs/gtdb_genomes_roary_eggnog/{gtdb_species}.emapper.annotations"
    conda: 'envs/eggnog.yml'
    resources:
        mem_mb = 32000,
        tmpdir=TMPDIR
    benchmark: "benchmarks/eggnog_{gtdb_species}.txt"
    threads: 8
    params: 
        outdir = "outputs/gtdb_genomes_roary_eggnog/",
        dbdir = "inputs/eggnog_db"
    shell:'''
    mkdir -p tmp/
    emapper.py --cpu {threads} -i {input.faa} --output {wildcards.gtdb_species} \
       --output_dir {params.outdir} -m hmmer -d none --tax_scope auto \
       --go_evidence non-electronic --target_orthologs all --seed_ortholog_evalue 0.001 \
       --seed_ortholog_score 60 --override --temp_dir tmp/ \
       -d 2 --data_dir {params.dbdir}
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

def checkpoint_charcoal_decontaminate_genomes3(wildcards):
    # checkpoint_output encodes the output dir from the checkpoint rule.
    checkpoint_output = checkpoints.charcoal_decontaminate_genomes.get(**wildcards).output[0]    
    file_names = expand("outputs/gtdb_genomes_intergenic_seqs/{{gtdb_species}}/{acc}.fa",
                        acc = glob_wildcards(os.path.join(checkpoint_output, "{acc}_genomic.fna.gz.clean.fa.gz")).acc)
    return file_names

rule cat_intergenic_sequences:
    input: checkpoint_charcoal_decontaminate_genomes3
    output: "outputs/gtdb_genomes_intergenic_comb/{gtdb_species}_all_intergenic_seqs.fa"
    threads: 1
    resources: 
        mem_mb=2000,
        tmpdir=TMPDIR
    benchmark:"benchmarks/cat_intergenic_{gtdb_species}.txt"
    shell:'''
    cat {input} > {output}
    '''

# note clustering might not be necessary because of the pufferfish index underlying salmon indexing. 
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
        seqs="outputs/gtdb_genomes_salmon_ref/{gtdb_species}.fa",
        decoys="outputs/gtdb_genomes_intergenic_comb/{gtdb_species}_clustered_intergenic_seq_names.txt"
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

            shell("fastp -i {params.tmp_base}.fastq.gz --json {output.json} --html {output.html} -R {wildcards.sra} --stdout > {output.reads}")

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

            shell("fastp -i {params.tmp_base}_1.fastq.gz -I {params.tmp_base}_2.fastq.gz --json {output.json} --html {output.html} -R {wildcards.sra} --stdout > {output.reads}")

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
        db="inputs/sourmash_dbs/gtdb-rs202.genomic.k31.zip",
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
        index = "outputs/gtdb_genomes_salmon_index/{gtdb_species}/info.json",
        reads = "outputs/rnaseq_fastp/{sra}.fq.gz"
    output: "outputs/rnaseq_salmon/{sra}/{gtdb_species}-{sra}_quant/quant.sf"
    params: 
        index_dir = lambda wildcards: "outputs/gtdb_genomes_salmon_index/" + wildcards.gtdb_species,
        out_dir = lambda wildcards: "outputs/rnaseq_salmon/{sra}/{gtdb_species_to_sra}_quant" 
    conda: "envs/salmon.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000 ,
        tmpdir= TMPDIR
    threads: 1
    benchmark: "benchmarks/rnaseq/salmon_quantify_{gtdb_species}-{sra}.txt"
    shell:'''
    salmon quant -i {params.index_dir} -l A -r {input.reads} -o {params.out_dir} --validateMappings
    '''
