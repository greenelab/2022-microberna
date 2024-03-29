import csv
import sys
import urllib.request
import pandas as pd

# set indir/outdir to read and write from
# on summit, indir is in projects, and outdir is scratch space
# since indir is the same directory that snakemake will be run from,
# it technically doesn't need to be parameterized
outdir = "/scratch/summit/treiter@xsede.org/2022-microberna/outputs"
indir = "/projects/treiter@xsede.org/2022-microberna/inputs"

m = pd.read_csv("inputs/20220407_runinfo.tsv.gz", sep = "\t", header = 0)
SRA = list(m['run_accession'])

gtdb_lineages = pd.read_csv("inputs/gtdb-rs207/gtdb-rs207.taxonomy.csv.gz", header = 0)
GTDB_ACCS = list(gtdb_lineages['ident'])

class Checkpoint_RnaseqToReference:
    """
    Define a class a la genome-grist to simplify file specification
    from checkpoint. This approach is documented at this url:
    http://ivory.idyll.org/blog/2021-snakemake-checkpoints.html
    """
    def __init__(self, pattern, sra_samples):
        self.pattern = pattern
        self.sra_samples = sra_samples


    def get_sra_to_species(self, sra):
        sra_to_ref_csv = f'{outdir}/rnaseq_sourmash_gather_to_ref_species/{sra}.csv'
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
        acc_csv = f'{outdir}/gtdb_genomes_by_species/{gtdb_species}.csv'
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
        # gather RNAseq sample
        #expand(outdir+"/rnaseq_sourmash_gather/{sra}_gtdb_k31.csv", sra = SRA),
        Checkpoint_RnaseqToReference(expand(outdir + "/rnaseq_salmon/{{gtdb_species}}/{sra}_quant/quant.sf", sra = SRA), SRA),

##############################################################
## Generate reference transcriptome using pangenome analysis
##############################################################

# The outputs of this rule are used by the class Checkpoint_GrabAccessions.
# The new wildcard, acc, is encoded in the "ident" column of the output
# file "accessions". 
checkpoint grab_species_accessions:
    input: gtdb_lineages = "inputs/gtdb-rs207/gtdb-rs207.taxonomy.csv.gz"
    output: accessions = outdir + "/gtdb_genomes_by_species/{gtdb_species}.csv"
    resources:
        mem_mb = 4000,
        time_min = 10 
    threads: 1
    conda: "envs/tidyverse.yml"
    script: "scripts/grab_species_accessions.R"


rule bakta_download_db:
    output: "inputs/bakta_db/db/version.json"
    threads: 1
    resources: 
        mem_mb = 4000,
        time_min = 240
    params: outdir = "inputs/bakta_db"
    conda: "envs/bakta.yml"
    shell:'''
    bakta_db download --output {params.outdir}
    '''

# TODO: add url for tar gz of genomes
rule download_charcoal_gtdb_rs207_genomes:
    output: "inputs/charcoal_gtdb_rs207_genomes.tar.gz"
    threads: 1
    resources: 
        mem_mb = 1000,
        time_min = 240
    shell:'''
    wget -O {output} URL
    '''

rule decompress_charcoal_gtdb_rs207_genomes:
    input: "inputs/charcoal_gtdb_rs207_genomes.tar.gz"
    output: expand("inputs/charcoal_gtdb_rs207_genomes/{gtdb_accs}_genomic.fna.gz.clean.fa.gz", gtdb_accs = GTDB_ACCS)
    params: outdir = "inputs/charcoal_gtdb_rs207_genomes/"
    resources: 
        mem_mb = 1000,
        time_min = 60
    shell:'''
    tar xf {input} -C {params.outdir}
    '''

rule bakta_annotate_gtdb_genomes:
    input: 
        fna=ancient("inputs/charcoal_gtdb_rs207_genomes/{acc}_genomic.fna.gz.clean.fa.gz"),
        db="inputs/bakta_db/db/version.json",
    output: 
        outdir + "/gtdb_genomes_bakta/{gtdb_species}/{acc}.faa",
        outdir + "/gtdb_genomes_bakta/{gtdb_species}/{acc}.gff3",
        outdir + "/gtdb_genomes_bakta/{gtdb_species}/{acc}.fna",
        outdir + "/gtdb_genomes_bakta/{gtdb_species}/{acc}.ffn",
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 16000 ,
        time_min = 40
    benchmark: "benchmarks/gtdb_genomes/bakta_{gtdb_species}/{acc}.txt"
    conda: 'envs/bakta.yml'
    params: 
        dbdir="inputs/bakta_db/db/",
        outdirb = lambda wildcards: f'{outdir}/gtdb_genomes_bakta/' + wildcards.gtdb_species,
    threads: 1
    shell:'''
    bakta --threads {threads} --db {params.dbdir} --prefix {wildcards.acc} --output {params.outdirb} \
        --locus {wildcards.acc} --locus-tag {wildcards.acc} --keep-contig-headers {input.fna}
    '''

rule cat_annotated_sequences:
    input: ancient(Checkpoint_GrabAccessions(outdir + "/gtdb_genomes_bakta/{{gtdb_species}}/{acc}.ffn"))
    output: outdir + "/gtdb_genomes_annotated_comb/{gtdb_species}_all_annotated_seqs.fa"
    threads: 1
    resources: 
        mem_mb=2000,
        time_min = 20
    benchmark:"benchmarks/gtdb_genomes/cat_annotated_seqs_{gtdb_species}.txt"
    shell:'''
    cat {input} > {output}
    '''

# note clustering might not be necessary because of the pufferfish index underlying salmon indexing. 
# but in any case it should help with downstream accounting, so keep.
rule cluster_annotated_sequences:
    input: outdir + "/gtdb_genomes_annotated_comb/{gtdb_species}_all_annotated_seqs.fa"
    output: outdir + "/gtdb_genomes_annotated_comb/{gtdb_species}_clustered_annotated_seqs.fa"
    threads: 1
    resources: 
        mem_mb=16000,
        time_min=30
    conda: "envs/cdhit.yml"
    benchmark:"benchmarks/gtdb_genomes/cluster_annotated_{gtdb_species}.txt"
    shell:'''
    cd-hit-est -c .95 -d 0 -i {input} -o {output}
    '''

rule translate_clustered_sequences_for_annotations:
    input: outdir + "/gtdb_genomes_annotated_comb/{gtdb_species}_clustered_annotated_seqs.fa" 
    output: outdir + '/gtdb_genomes_annotated_comb/{gtdb_species}_clustered_annotated_seqs.faa'
    conda: 'envs/emboss.yml'
    resources:
        mem_mb = 16000,
        time_min=20
    benchmark: "benchmarks/gtdb_genomes/transeq_{gtdb_species}.txt"
    threads: 2
    shell:'''
    transeq {input} {output}
    '''

rule eggnog_download_db:
    output: "inputs/eggnog_db/eggnog.db"
    threads: 1   
    resources: 
        mem_mb = 4000,
        time_min=2880
    params: datadir = "inputs/eggnog_db"
    conda: "envs/eggnog.yml"
    shell:'''
    download_eggnog_data.py -H -d 2 -y --data_dir {params.datadir}
    '''

rule eggnog_annotate_clustered_sequences:
    input: 
        faa = outdir + "/gtdb_genomes_annotated_comb/{gtdb_species}_clustered_annotated_seqs.faa",
        db = 'inputs/eggnog_db/eggnog.db'
    output: outdir + "/gtdb_genomes_annotated_comb_eggnog/{gtdb_species}/{gtdb_species}.emapper.annotations"
    conda: 'envs/eggnog.yml'
    resources:
        mem_mb = 32000,
        time_min=2880
    benchmark: "benchmarks/gtdb_genomes/eggnog_{gtdb_species}.txt"
    threads: 8
    params: 
        outdire = lambda wildcards: outdir + "/gtdb_genomes_annotated_comb_eggnog/" + wildcards.gtdb_species,
        dbdir = "inputs/eggnog_db"
    shell:'''
    mkdir -p tmp/
    emapper.py --cpu {threads} -i {input.faa} --output {wildcards.gtdb_species} \
       --output_dir {params.outdire} -m hmmer -d none --tax_scope auto \
       --go_evidence non-electronic --target_orthologs all --seed_ortholog_evalue 0.001 \
       --seed_ortholog_score 60 --override --temp_dir tmp/ \
       -d 2 --data_dir {params.dbdir}
    '''

#########################################################
## generate pangenome
#########################################################

rule roary_determine_pangenome:
    input: ancient(Checkpoint_GrabAccessions(outdir + "/gtdb_genomes_bakta/{{gtdb_species}}/{acc}.gff3"))
    output: outdir + "/gtdb_genomes_roary/{gtdb_species}/pan_genome_reference.fa"
    conda: 'envs/roary.yml'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000,
        time_min = 2880
    benchmark: "benchmarks/gtdb_genomes/roary_{gtdb_species}.txt"
    threads: 2
    params: 
        outdirr = lambda wildcards: outdir + '/gtdb_genomes_roary/' + wildcards.gtdb_species
    shell:'''
    roary -r -e -n -f {params.outdirr} -p {threads} -z {input}
    '''


############################################################
## Generate decoy sequences
############################################################

rule define_genome_regions_as_bed:
    input: ancient(outdir + "/gtdb_genomes_bakta/{gtdb_species}/{acc}.fna"),
    output: outdir + "/gtdb_genomes_regions_bed/{gtdb_species}/{acc}_region.sizes"
    threads: 1
    resources: 
        mem_mb=2000,
        tim_min = 10
    conda: "envs/samtools.yml"
    benchmark:"benchmarks/gtdb_genomes/genome_regions_{gtdb_species}/{acc}.txt"
    shell:'''    
    samtools faidx {input}
    cut -f 1,2 {input}.fai | sort -n -k1 > {output}
    '''

rule filter_bakta_gff:
    input:  gff=ancient(outdir + "/gtdb_genomes_bakta/{gtdb_species}/{acc}.gff3"),
    output: gff=temp(outdir + "gtdb_genomes_gff/{gtdb_species}/{acc}_filtered.gff3"),
    threads: 1
    resources: 
        mem_mb=2000,
        time_min = 10
    conda: "envs/rtracklayer.yml"
    benchmark:"benchmarks/gtdb_genomes/filter_gff_{gtdb_species}/{acc}.txt"
    script: "scripts/filter_gff.R"

rule sort_bakta_gff:
    input: gff=outdir + "/gtdb_genomes_gff/{gtdb_species}/{acc}_filtered.gff3",
    output: gff=outdir + "/gtdb_genomes_gff/{gtdb_species}/{acc}_filtered_sorted.gff3",
    threads: 1
    resources: 
        mem_mb=2000,
        time_min=5
    conda: "envs/bedtools.yml"
    benchmark:"benchmarks/gtdb_genomes/sort_gff_{gtdb_species}/{acc}.txt"
    shell:'''
    sortBed -i {input} > {output}
    '''

rule complement_gff:
    input:
        gff=outdir + "/gtdb_genomes_gff/{gtdb_species}/{acc}_filtered_sorted.gff3",
        sizes= outdir + "/gtdb_genomes_regions_bed/{gtdb_species}/{acc}_region.sizes"
    output: outdir + "/gtdb_genomes_intergenic_bed/{gtdb_species}/{acc}.bed"
    threads: 1
    resources: 
        mem_mb=2000,
        time_min=5
    conda: "envs/bedtools.yml"
    benchmark:"benchmarks/gtdb_genomes/complement_gff_{gtdb_species}/{acc}.txt"
    shell:'''
    bedtools complement -i {input.gff} -g {input.sizes} > {output}
    '''

rule extract_intergenic_from_genome:
    input:
        fna=ancient(outdir+"/gtdb_genomes_bakta/{gtdb_species}/{acc}.fna"),
        bed=outdir+"/gtdb_genomes_intergenic_bed/{gtdb_species}/{acc}.bed"
    output: outdir + "/gtdb_genomes_intergenic_seqs/{gtdb_species}/{acc}.fa"
    threads: 1
    resources: 
        mem_mb=2000,
        time_min=5
    conda: "envs/bedtools.yml"
    benchmark:"benchmarks/gtdb_genomes/extract_intergenic_{gtdb_species}/{acc}.txt"
    shell:'''
    bedtools getfasta -fi {input.fna} -bed {input.bed} -name -s > {output}
    '''

rule cat_intergenic_sequences:
    input: ancient(Checkpoint_GrabAccessions(outdir+"/gtdb_genomes_intergenic_seqs/{{gtdb_species}}/{acc}.fa")),
    output: outdir+"/gtdb_genomes_intergenic_comb/{gtdb_species}_all_intergenic_seqs.fa"
    threads: 1
    resources: 
        mem_mb=2000,
        time_min=5
    benchmark:"benchmarks/gtdb_genomes/cat_intergenic_{gtdb_species}.txt"
    shell:'''
    cat {input} > {output}
    '''

rule cluster_intergenic_sequences:
    input: outdir+"/gtdb_genomes_intergenic_comb/{gtdb_species}_all_intergenic_seqs.fa"
    output: outdir+"/gtdb_genomes_intergenic_comb/{gtdb_species}_clustered_intergenic_seqs.fa"
    threads: 1
    resources: 
        mem_mb=16000,
        time_min=30
    conda: "envs/cdhit.yml"
    benchmark:"benchmarks/gtdb_genomes/cluster_intergenic_{gtdb_species}.txt"
    shell:'''
    cd-hit-est -c 1 -i {input} -o {output}
    '''

###################################################
## Generate decoy-aware salmon transcriptome index
###################################################

rule grab_intergenic_sequence_names:
    input: outdir+"/gtdb_genomes_intergenic_comb/{gtdb_species}_clustered_intergenic_seqs.fa"
    output: outdir+"/gtdb_genomes_intergenic_comb/{gtdb_species}_clustered_intergenic_seq_names.txt"
    threads: 1
    resources: 
        mem_mb=2000,
        time_min=10
    benchmark:"benchmarks/gtdb_genomes/grab_intergenic_names_{gtdb_species}.txt"
    shell:"""
    grep '>' {input} | awk 'sub(/^>/, "")' > {output}
    """

rule combine_pangenome_and_intergenic_sequences:
    input: 
        outdir+"/gtdb_genomes_annotated_comb/{gtdb_species}_clustered_annotated_seqs.fa",
        outdir+"/gtdb_genomes_intergenic_comb/{gtdb_species}_clustered_intergenic_seqs.fa"
    output: outdir+"/gtdb_genomes_salmon_ref/{gtdb_species}.fa"
    threads: 1
    resources: 
        mem_mb=2000,
        time_min=5
    benchmark:"benchmarks/gtdb_genomes/combine_sequences_{gtdb_species}.txt"
    shell:'''
    cat {input} > {output}
    '''

rule index_transcriptome:
    input: 
        seqs=ancient(outdir+"/gtdb_genomes_salmon_ref/{gtdb_species}.fa"),
        decoys=ancient(outdir+"/gtdb_genomes_intergenic_comb/{gtdb_species}_clustered_intergenic_seq_names.txt")
    output: outdir+"/gtdb_genomes_salmon_index/{gtdb_species}/info.json"
    threads: 1
    params: index_dir = lambda wildcards: outdir+"/gtdb_genomes_salmon_index/" + wildcards.gtdb_species
    resources: 
        mem_mb=16000,
        time_min=20
    conda: "envs/salmon.yml"
    benchmark:"benchmarks/gtdb_genomes/salmon_index_{gtdb_species}.txt"
    shell:'''
    salmon index -t {input.seqs} -i {params.index_dir} --decoys {input.decoys} -k 31
    '''

#################################################################
## Download, process, and taxonomically annotate RNAseq samples
#################################################################

rule rnaseq_sample_download:
    output:
        reads=outdir+"/rnaseq_fastp/{sra}.fq.gz",
        json = outdir+"/rnaseq_fastp/{sra}.fastp.json",
        html = outdir+"/rnaseq_fastp/{sra}.fastp.html"
    params: tmp_base = lambda wildcards: outdir + "/tmp_raw/" + wildcards.sra
    threads: 1
    resources:
        mem_mb=8000,
        time_min=960
    run:
        row = m.loc[m['run_accession'] == wildcards.sra]
        fastqs = row['fastq_ftp'].values[0]
        fastqs = fastqs.split(";")
        if len(fastqs) == 1:
            # single end data; download and stream directly to fastp for trimming.
            fastq = fastqs[0]
            shell("mkdir -p {outdir}/tmp_raw")
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
            shell("mkdir -p {outdir}/tmp_raw")
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
    input: outdir+"/rnaseq_fastp/{sra}.fq.gz",
    output: outdir+"/rnaseq_sourmash_sketch/{sra}.sig"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000 ,
        time_min=120
    threads: 1
    benchmark: "benchmarks/rnaseq/sourmash_sketch_{sra}.txt"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sketch dna -p k=21,k=31,k=51,scaled=1000,abund -o {output} --name {wildcards.sra} {input}
    '''

# TODO: update download link to RS207 when on OSF
rule sourmash_download_gtdb_database:
    output: "inputs/sourmash_dbs/gtdb-rs207.genomic.dna.k31.zip"
    resources:
        mem_mb = 1000,
        time_min = 240
    threads: 1
    shell:'''
    wget -O {output} https://osf.io/94mzh/download
    '''

rule sourmash_download_human_sig:
    output: "inputs/sourmash_dbs/GCF_000001405.39_GRCh38.p13_rna.sig"
    resources:
        mem_mb = 1000,
        time_min = 10
    threads: 1
    shell:'''
    wget -O {output} https://osf.io/anj6b/download
    '''

rule sourmash_download_euk_database:
    output: "inputs/sourmash_dbs/euk_rna_k31.tar.gz"
    resources:
        mem_mb = 1000,
        time_min = 10
    threads: 1
    shell:'''
    wget -O {output} https://osf.io/qk5th/download
    '''

rule sourmash_euk_database_untar:
    '''
    euk rna is a legacy database. We've since switched
    formats, but this one needs to be untarred before it 
    can be used for gather
    '''
    input: "inputs/sourmash_dbs/euk_rna_k31.tar.gz"
    output: "inputs/sourmash_dbs/euk_rna_k31.sbt.json"
    params: out_dir = "inputs/sourmash_dbs"
    resources:
        mem_mb = 1000,
        time_min = 10
    threads: 1
    shell:'''
    tar xf {input} -C {params.out_dir}
    '''

rule rnaseq_sample_sourmash_gather_against_gtdb:
    input:
        sig=outdir+"/rnaseq_sourmash_sketch/{sra}.sig",
        db="inputs/sourmash_dbs/gtdb-rs207.genomic.dna.k31.zip",
        human="inputs/sourmash_dbs/GCF_000001405.39_GRCh38.p13_rna.sig",
        euk="inputs/sourmash_dbs/euk_rna_k31.sbt.json"
    output: outdir+"/rnaseq_sourmash_gather/{sra}_gtdb_k31.csv"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 34000 ,
        time_min = 180
    threads: 1
    benchmark: "benchmarks/rnaseq/sourmash_gather_k31_{sra}.txt"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash gather -o {output} --scaled 2000 -k 31 {input.sig} {input.db} {input.human} {input.euk}
    '''

rule fix_dag:
    input: expand(outdir+"/rnaseq_sourmash_gather/{sra}_gtdb_k31.csv", sra =SRA)
    output: touch(outdir + "/rnaseq_sourmash_gather/all_gather_done.txt")

checkpoint rnaseq_sample_select_best_species_reference:
    """
    Corresponds to class Checkpoint_RnaseqToReference at top of file.
    This will read in the gather results, and output one dataframe per RNA seq sample. 
    The dataframe will have the SRA accession and the reference genome accession.
    wildcard {sra} will already be defined as a list of RNAseq SRA accessions at the
    beginning of the Snakefile, and the checkpoint class will output the proper 
    gtdb species for each SRA accession.
    """
    input:
        dummy = outdir + "/rnaseq_sourmash_gather/all_gather_done.txt",
        gtdb_lineages="inputs/gtdb-rs207/gtdb-rs207.taxonomy.csv.gz",
        gather=ancient(outdir+"/rnaseq_sourmash_gather/{sra}_gtdb_k31.csv")
    output: sra_to_ref_species= outdir + "/rnaseq_sourmash_gather_to_ref_species/{sra}.csv"
    conda: "envs/tidyverse.yml"
    resources:
        mem_mb = 2000,
        time_min=5
    threads: 1
    benchmark: "benchmarks/rnaseq/select_species_genome_{sra}.txt"
    script: "scripts/select_best_species_reference.R"

rule rnaseq_quantify_against_species_pangenome:
    input: 
        sra_to_ref_species=ancient(outdir+"/rnaseq_sourmash_gather_to_ref_species/{sra}.csv"),
        index = ancient(outdir+"/gtdb_genomes_salmon_index/{gtdb_species}/info.json"),
        reads = outdir+"/rnaseq_fastp/{sra}.fq.gz",
        eggnog = outdir+"/gtdb_genomes_annotated_comb_eggnog/{gtdb_species}/{gtdb_species}.emapper.annotations" # not actually needed here, but trick snakemake into making these annots
    output: outdir+"/rnaseq_salmon/{gtdb_species}/{sra}_quant/quant.sf"
    params: 
        index_dir = lambda wildcards: outdir + "/gtdb_genomes_salmon_index/" + wildcards.gtdb_species,
        out_dir = lambda wildcards: outdir + "/rnaseq_salmon/" + wildcards.gtdb_species + "/" + wildcards.sra + "_quant" 
    conda: "envs/salmon.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000 ,
        time_min=180
    threads: 1
    benchmark: "benchmarks/rnaseq/salmon_quantify_{gtdb_species}-{sra}.txt"
    shell:'''
    salmon quant -i {params.index_dir} -l A -r {input.reads} -o {params.out_dir} --validateMappings --writeUnmappedNames
    '''
