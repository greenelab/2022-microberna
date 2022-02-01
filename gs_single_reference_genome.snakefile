# DYNAMIC APPROACH
# download and qc reads
# run gather
# for samples that have a single gather reference, select that reference, download the genome, and index


# NON-DYNAMIC APPROACH
# define list of known single-reference samples (e.g. Fpc -- all are dominantly the same strain)
# use alpha ksize approach to bind and limit wildcards (can probably use pandas to make these lists from a csv file that records SRA number and reference genome ID)
# download reference, reference gff file, and index with star
# download and qc reads
# map with star
# count with e.g. fadu

import pandas as pd
import re
import sys
import csv

TMPDIR = "/tmp"

m = pd.read_csv("inputs/toy_metadata/gs_fpc.tsv", header = 0, sep =  "\t")
m["genome_acc_to_sra"] = m["reference_genome_accession"] + "-" + m["experiment_accession"]
SRA = m["experiment_accession"]
GENOME_ACC = m["reference_genome_accession"]
# bind each SRA accession to its reference genome accession to limit combinations of wildcards
GENOME_ACC_TO_SRA = m["genome_acc_to_sra"]

rule all:
    input:
        expand('outputs/gs_rnaseq_star/{acc_to_sra}Aligned.sortedByCoord.out.bam', acc_to_sra = GENOME_ACC_TO_SRA) 

rule make_genome_info_csv:
    output:
        csv = 'inputs/gs_genomes/{acc}.info.csv'
    conda: "envs/download_genomes.yml"
    resources:
        mem_mb = 8000
    threads: 1
    shell: """
    python scripts/genbank_genomes.py {wildcards.acc} --output {output.csv}
    """

rule download_matching_genomes:
    input:
        csvfile = ancient('inputs/gs_genomes/{acc}.info.csv')
    output:
        genome = "inputs/gs_genomes/{acc}_genomic.fna.gz"
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

rule bakta_download_db:
    output: "inputs/bakta_db/db/version.json"
    threads: 1
    resources: mem_mb = 4000
    params: outdir = "inputs/bakta_db"
    conda: "envs/bakta.yml"
    shell:'''
    bakta_db download --output {params.outdir}
    '''

rule bakta_annotate_gs_genomes:
    # TODO: change locus_tag to accept wildcards.acc if bakta/#92 gets implemented.
    input:
        fna=ancient("inputs/gs_genomes/{acc}_genomic.fna.gz"),
        db="inputs/bakta_db/db/version.json",
    output:
        "outputs/gs_genomes_bakta/{acc}/{acc}.faa",
        "outputs/gs_genomes_bakta/{acc}/{acc}.gff3",
        "outputs/gs_genomes_bakta/{acc}/{acc}.fna",
        "outputs/gs_genomes_bakta/{acc}/{acc}.ffn",
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000 ,
        tmpdir= TMPDIR
    benchmark: "benchmarks/gs_single_reference_genome/bakta_{acc}.txt"
    conda: 'envs/bakta.yml'
    params:
        dbdir="inputs/bakta_db/db/",
        outdir = lambda wildcards: 'outputs/gs_genomes_bakta/' + wildcards.acc,
        #locus_tag = lambda wildcards: re.sub("[CF_\.]", "", wildcards.acc)
    threads: 1
    shell:'''
    bakta --threads {threads} --db {params.dbdir} --prefix {wildcards.acc} --output {params.outdir} \
        --locus {wildcards.acc} --locus-tag {wildcards.acc} --keep-contig-headers {input.fna}
    '''

rule convert_gff_to_gtf_gs_genomes:
    input: gff="outputs/gs_genomes_bakta/{acc}/{acc}.gff3",
    output: gtf="outputs/gs_genomes_bakta/{acc}/{acc}.gtf",
    conda: "envs/rtracklayer.yml"
    threads: 1
    resources: 
        mem_mb = 3000,
        tmpdir= TMPDIR
    script: "scripts/convert_gff_to_gtf.R"

rule rnaseq_sample_download:
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

rule star_index_genome:
    input:
        fna = "outputs/gs_genomes_bakta/{acc}/{acc}.fna",
    output: 'outputs/gs_genomes_bakta/{acc}/SAindex'
    params: input_dir = 'outputs/gs_genomes_bakta/{acc}/' 
    conda: 'envs/star.yml'
    resources: 
        mem_mb = 16000,
        tmpdir= TMPDIR
    shell:'''
    STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {params.input_dir} \
         --genomeFastaFiles {input.fna} --genomeSAindexNbases 8
    '''

rule star_map_single_end:
    input:
        reads = "outputs/rnaseq_fastp/{sra}.fq.gz",
        genome_index = 'outputs/gs_genomes_bakta/{acc}/SAindex'
    output: 'outputs/gs_rnaseq_star/{acc}-{sra}Aligned.sortedByCoord.out.bam' 
    params: 
        out_prefix = lambda wildcards: 'outputs/gs_rnaseq_star/' + wildcards.acc + "-" + wildcards.sra,
        genome_dir = lambda wildcards: 'outputs/gs_genomes_bakta/' + wildcards.acc
    conda: 'envs/star.yml' 
    threads: 2   
    resources: 
        mem_mb = 16000,
        tmpdir= TMPDIR
    shell:'''
    STAR --runThreadN {threads} --genomeDir {params.genome_dir}      \
        --readFilesIn {input.reads} --outFilterMultimapNmax 20 \
        --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverLmax 0.6 --alignIntronMin 20 \
        --outSAMattributes NH HI NM MD --outSAMtype BAM SortedByCoordinate \
        --limitBAMsortRAM 1040169341 --outFileNamePrefix {params.out_prefix}
    '''

rule count_single_end:
    input: 
        bam='outputs/gs_rnaseq_star/{acc}-{sra}Aligned.sortedByCoord.out.bam',
        gtf="outputs/gs_genomes_bakta/{acc}/{acc}.gtf",
    output:
        counts="outputs/gs_rnaseq_featurecounts/{acc}-{sra}_counts.txt",
        juncs="outputs/gs_rnaseq_featurecounts/{acc}-{sra}_junctions.txt"
    conda: 'envs/rsubread.yml' 
    resources: 
        mem_mb = 8000,
        tmpdir= TMPDIR
    threads: 1
    script: "scripts/featureCounts_se.R"
      

# not all sequences are paired end; decide later if pe should be benchmarked explicitly
#rule split_paired_end_reads:
#rule star_map_paired_end:
#rule count_paired_end:
