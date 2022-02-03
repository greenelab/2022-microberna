import csv
import sys
import urllib.request
import pandas as pd

m = pd.read_csv("inputs/toy_metadata/toy_fpc.tsv", sep = "\t", header = 0)
SRA = list(m['experiment_accession'])

TMPDIR = "/scratch/tereiter/" # TODO: update tmpdir based on computing env, or remove tmpdir invocation in resources

class Checkpoint_RnaseqToReference:
    """
    Define a class a la genome-grist to simplify file specification
    from checkpoint (e.g. solve for {gtdb_species} wildcard). This approach
    is documented at this url:
    http://ivory.idyll.org/blog/2021-snakemake-checkpoints.html
    """
    def __init__(self, pattern, sra_samples):
        self.pattern = pattern
        self.sra_samples = sra_samples
        print('XXX', pattern)
        print('YYY', (sra_samples,), type(sra_samples))

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

        print('ZZZ', results)

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
        #expand("outputs/rnaseq_sourmash_gather/{sra}_gtdb_k31.csv", sra = SRA)
        Checkpoint_RnaseqToReference(expand("outputs/rnaseq_salmon/{{gtdb_species}}/{sra}_quant/quant.sf", sra = SRA), SRA)

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


rule cat_genome_sequences:
    input: Checkpoint_GrabAccessions("inputs/gtdb_genomes/{{gtdb_species}}/{acc}_genomic.fna.gz")
    output: "outputs/gtdb_genomes_cat/{gtdb_species}.fna.gz"
    threads: 1
    resources: 
        mem_mb=2000,
        tmpdir=TMPDIR
    benchmark:"benchmarks/cat_annotated_seqs_{gtdb_species}.txt"
    shell:'''
    cat {input} > {output}
    '''

rule index_genome_sequences:
    input: 
        seqs="outputs/gtdb_genomes_cat/{gtdb_species}.fna.gz",
    output: "outputs/gtdb_genomes_salmon_index/{gtdb_species}/info.json"
    threads: 1
    params: index_dir = lambda wildcards: "outputs/gtdb_genomes_salmon_index/" + wildcards.gtdb_species
    resources: 
        mem_mb=16000,
        tmpdir=TMPDIR
    conda: "envs/salmon.yml"
    benchmark:"benchmarks/salmon_index_{gtdb_species}.txt"
    shell:'''
    salmon index -t {input.seqs} -i {params.index_dir} -k 31
    '''

#################################################################
## Download, process, and taxonomically annotate RNAseq samples
#################################################################

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
                shell("curl ftp://{fastq} | zcat | head -n 1000 | gzip > {params.tmp_base}.fastq.gz")

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
                shell("curl ftp://{fastq_1} | zcat | head -n 500 | gzip > {params.tmp_base}_1.fastq.gz")

            if not os.path.exists(params.tmp_base + "_2.fastq.gz"):
                shell("curl ftp://{fastq_2} | zcat | head -n 500 | gzip > {params.tmp_base}_2.fastq.gz")

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

rule sourmash_download_gtdb_database:
    output: "inputs/sourmash_dbs/gtdb-rs202.genomic.k31.zip" 
    resources:
        mem_mb = 1000,
        tmpdir = TMPDIR
    threads: 1
    shell:'''
    wget -O {output} https://osf.io/94mzh/download
    '''

rule rnaseq_sample_sourmash_gather_against_gtdb:
    input:
        sig="outputs/rnaseq_sourmash_sketch/{sra}.sig",
        db="inputs/sourmash_dbs/gtdb-rs202.genomic.k31.zip",
    output: "outputs/rnaseq_sourmash_gather/{sra}_gtdb_k31.csv"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000 ,
        tmpdir= TMPDIR
    threads: 1
    benchmark: "benchmarks/rnaseq/sourmash_gather_k31_{sra}.txt"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash gather -o {output} --scaled 2000 -k 31 {input.sig} {input.db}
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
    salmon quant -i {params.index_dir} -l A -r {input.reads} -o {params.out_dir} --validateMappings
    '''
