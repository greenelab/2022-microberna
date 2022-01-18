import gzip
import sys
from os import path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
#import pyrodigal
import pandas as pd

TMPDIR = "/tmp"

m = pd.read_csv("inputs/toy_metadata/toy_pa.tsv", sep = "\t", header = 0)
SRA = m.sort_values(by='read_count')['experiment_accession']

# tmp PA samples:
# PAO1: ERX3558803, SRX540112
# clinical isolates: SRX396878, labelled B271; SRX5123759
# PA14: SRX589549, SRX4624095

# TODO: update file/checkpoint names
# TODO: switch out column names
# TODO: add in wildcard for SRA_accession (experiment or sample accession? experiment has more in my SRA df)
class Checkpoint_RnaseqToReference:
    """
    Define a class a la genome-grist to simplify file specification
    from checkpoint (e.g. solve for {acc} wildcard). This approach
    is documented at this url:
    http://ivory.idyll.org/blog/2021-snakemake-checkpoints.html
    """
    def __init__(self, pattern):
        self.pattern = pattern

    def get_acc_dbs(self):
        acc_db_csv = f'outputs/genbank/gather_gtdb-rs202-genomic.x.dbs.csv'
        assert os.path.exists(acc_db_csv)

        acc_dbs = []
        with open(acc_db_csv, 'rt') as fp:
           r = csv.DictReader(fp)
           for row in r:
               acc = row['accession']
               db = row['species']
               acc_db = acc + "-" + db
               acc_dbs.append(acc_db)

        return acc_dbs

    def __call__(self, w):
        global checkpoints

        # wait for the results of 'query_to_species_db';
        # this will trigger exception until that rule has been run.
        checkpoints.query_to_species_db.get(**w)

        # parse accessions in gather output file
        genome_acc_dbs = self.get_acc_dbs()

        p = expand(self.pattern, acc_db=genome_acc_dbs, **w)
        return p


rule all:
    input:
        expand("outputs/rnaseq_sourmash_gather/{sra}_gtdb_k31.csv", sra = SRA)

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
        mem_mb=1000
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
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sketch dna -p k=21,k=31,scaled=2000 -o {output} --name {wildcards.sra} {input}
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
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash gather -o {output} --scaled 2000 -k 31 {input.sig} {input.db} {input.human}
    '''

checkpoint rnaseq_sample_select_best_genome_reference:
    """
    I think this will read in the gather results, and output 
    one dataframe per RNA seq sample. The dataframe will have the SRA
    accession and the reference genome accession. wildcard {sra} will already
    be defined as a list of RNAseq SRA accessions at the beginning of the 
    Snakefile, and the checkpoint class will output an amalgamation of wildcards
    {sra}-{acc}, which will be called {sra_to_acc} in the rule all, but the hyphen
    should allow for the wildcards to be solved separately. 
    Example:
    https://github.com/taylorreiter/2021-metapangenome-example/blob/main/Snakefile#L90
    """

rule reference_make_accession_info_csv:
    """
    Given an accession, find the NCBI GenBank url to download the associated data products.
    """
    output:
        csvfile = 'inputs/gtdb_genomes/{acc}.info.csv'
    shell: """
    python -Werror scripts/genbank_cds.py {wildcards.acc} \
       --output {output.csvfile}
    """

rule reference_download_matching_cds_or_download_genome_and_predict_cds:
    """
    Download cds_from_genomic.fna.gz file if it exists.  
    If not, download the genome and use pyrodigal (python wrapper for
    prodigal) to predict open reading frames in the genome.
    
    This rule was inspired by N. T. Pierce-Wards genome download rules in this snakefile:
    https://github.com/bluegenes/2021-rank-compare/blob/e463577851e6fbccb4a06d04fedaa430fedbe804/picklist-bf.snakefile
    The rules in the linked snakefile were in turn inspired by the 
    genome download rules by C. Titus Brown in the genome-grist pipeline. 
    """
    input:
        csvfile = ancient('genbank/info/{acc}.info.csv')
    output:
        cds = "outputs/gtdb_genome_cds/{acc}_cds_from_genomic.fna.gz"
    params: fna = lambda wildcards: "outputs/gtdb_genome_fna/" + wildcards.acc + "_genomic.fna.gz"
    run:
        with open(input.csvfile, 'rt') as infp:
            r = csv.DictReader(infp)
            rows = list(r)
            assert len(rows) == 1
            row = rows[0]
            acc = row['acc']
            assert wildcards.acc.startswith(acc)
            cds_url = row['cds_url']
            genome_url = row['genome_url']
            name = row['ncbi_tax_name']

            print(f"downloading cds for acc {acc}/{name} from GenBank...",
                file=sys.stderr)
            try:
                with open(output.cds, 'wb') as outfp:
                    with urllib.request.urlopen(cds_url) as response:
                        content = response.read()
                        outfp.write(content)
                        print(f"...wrote {len(content)} bytes to {output.cds}",
                              file=sys.stderr)
            except:
                # download the genome if the cds file doesn't exist
                with open(params.fna, 'wb') as outfp:
                    with urllib.request.urlopen(genome_url) as response:
                        content = response.read()
                        outfp.write(content)
                        print(f"...wrote {len(content)} bytes to genome file")
                # use pyrodigal to read in genomic file and predict cds from it
                with gzip.open(params.fna, "rt") as handle:
                    record = SeqIO.read(handle, 'fasta') # read in genome fasta file
                
                p = pyrodigal.OrfFinder(meta=True)
                predicted_orfs = []
                for prediction in p.find_genes(bytes(record.seq)):
                    predicted_orf = prediction.sequence()
                    predicted_orfs.append(predicted_orf)

                with gzip.open(output.cds, "wt") as file:
                    print("writing pyrodigal predicted cds...")
                    for (index, seq) in enumerate(predicted_orfs):
                    # prepare row
                        row = SeqRecord(Seq(seq), id = record.name + "_" + str(index), name = record.name, description = "")
                        SeqIO.write(row, file, "fasta")
                    


# rule tmp unzip cds
# rule emboss transeq
# rule eggnog annotate emboss transeq
# rule reference_salmon_index:
# rule rna_seq_sample_salmon_quant:
#
# IDEAS
# + compendia of all rna-seq that are the same species using the NCBI gene identifiers (or whatever refinebio is already using for e.g. e. coli k12)
# + species-level compendia annotated by pangenome? or annotated by ortholog identifiers
#     + pangenome would be download all genomes, annotate with bakta, run roary, assign identifiers and collect
#     + annotate by ortholog identifers would be CDS -> aa seqs -> eggnog 
