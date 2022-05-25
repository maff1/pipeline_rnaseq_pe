#########################################
# Snakemake pipeline for RNA-Seq analysis
#########################################

###############
# libraries
###############
import os.path
import pandas as pd

###############
# Configuration
###############
configfile: "config/config.yaml"
REFERENCE_DIR = config["reference_dir"]
RESULT_DIR = config["result_dir"]

samples_table = pd.read_csv(config["samples_tbl"]).set_index("sample", drop=False)
SAMPLES = samples_table.index.get_level_values('sample').unique().tolist()

###########################
# Input functions for rules
###########################

def fq1_sample(wildcards):
    return samples_table.loc[wildcards.sample, "fq1"]
    
def fq2_sample(wildcards):
    return samples_table.loc[wildcards.sample, "fq2"]

###########################
# Rules
###########################

rule all:
    input:
        REFERENCE_DIR + os.path.basename(config["refs"]["ensemblGTF"]),
        REFERENCE_DIR + os.path.basename(config["refs"]["ensemblFASTA"]),
        REFERENCE_DIR + os.path.basename(config["refs"]["ensembl_transc_fasta"]),
        REFERENCE_DIR + "Homo_sapiens.GRCh38.cdna.all.fa.gz",
        RESULT_DIR + "multiqc/multiqc_report.htm",
        REFERENCE_DIR + "salmon_index",
        REFERENCE_DIR + "genome/",
        RESULT_DIR + "mapping_summary.csv",
        RESULT_DIR + "tximeta/" + "tximeta.rds",
        RESULT_DIR + "gene_expression/" + "dds.rds",
        RESULT_DIR + "gene_expression/" + "dge.xlsx",
        RESULT_DIR + "gene_expression/" + "counts.tsv"
        
        
        

rule download_ensembl:
    params:
        gtf = config["refs"]["ensemblGTF"],
        fasta = config["refs"]["ensemblFASTA"],
        tfa = config["refs"]["ensembl_transc_fasta"]
    output:
        gtf_out = REFERENCE_DIR + os.path.basename(config["refs"]["ensemblGTF"]),
        fasta_out = REFERENCE_DIR + os.path.basename(config["refs"]["ensemblFASTA"]),
        tfa_out = REFERENCE_DIR + os.path.basename(config["refs"]["ensembl_transc_fasta"])
    shell:
        """
        wget --no-check-certificate {params.gtf} -O {output.gtf_out}; \
        wget --no-check-certificate {params.fasta} -O {output.fasta_out}; \
        wget --no-check-certificate {params.tfa} -O {output.tfa_out};
        """

rule salmon_indexing:
    input:
        REFERENCE_DIR + "Homo_sapiens.GRCh38.cdna.all.fa.gz"
    threads: 12
    conda: "envs/salmon.yaml"
    output:
        directory(REFERENCE_DIR + "salmon_index")
    shell:
        """
        salmon index -t <(gunzip -c {input}) -p {threads} -i {output} -k 31
        """

#######################
# RNA-seq read trimming
#######################

rule qc_fastp:
    input:
        fq1 = fq1_sample,
        fq2 = fq2_sample
    output:
        fq1  = RESULT_DIR + "trimmed/" + "{sample}_R1_trimmed.fq.gz", 
        fq2  = RESULT_DIR + "trimmed/" + "{sample}_R2_trimmed.fq.gz",
        html = RESULT_DIR + "fastp/{sample}_fastp.html",
        json = RESULT_DIR + "fastp/{sample}_fastp.json",
    message:"trimming {wildcards.sample} reads"
    threads: 10
    log:
        RESULT_DIR + "fastp/{sample}.log.txt"
    resources: cpus=10
    conda: "envs/fastp.yaml"
    shell:
        "touch {output.fq2};\
        fastp --in1 {input.fq1} --in2 {input.fq2} --out1 {output.fq1} --out2 {output.fq2} \
        --thread {threads} --dedup --html {output.html} --json {output.json} \
        2>{log}"

rule multiqc:
    input:
        fastp = expand(RESULT_DIR + "fastp/{sample}_fastp.json", sample=SAMPLES),
        star = expand(RESULT_DIR + "star/{sample}_Log.final.out", sample=SAMPLES)
    output:
        directory(RESULT_DIR + "multiqc/multiqc_report.htm")
    conda: "envs/multiqc.yaml"
    shell:
        "multiqc {input.fastp} {input.star} -o {output}"

rule salmon_quant:
    input:
        fq1  = expand(RESULT_DIR + "trimmed/" + "{sample}_R1_trimmed.fq.gz", sample=SAMPLES),
        fq2  = expand(RESULT_DIR + "trimmed/" + "{sample}_R2_trimmed.fq.gz", sample=SAMPLES)
    output:
        directory(RESULT_DIR + "salmon_quant/" + "{sample}")
    params:
        index = REFERENCE_DIR + "salmon_index"
    threads: 12
    conda: "envs/salmon.yaml"
    shell:
        """
        salmon quant -i {params.index} -l A \
        -1 {input.fq1} \
        -2 {input.fq2} \
        -p {threads} --validateMappings --seqBias --gcBias --posBias -o {output}
        """
        
rule star_index:
    input:
        fasta = REFERENCE_DIR + os.path.basename(config["refs"]["ensemblFASTA"]),
        gtf = REFERENCE_DIR + os.path.basename(config["refs"]["ensemblGTF"])
    output: directory(REFERENCE_DIR + "genome/")
    params:
        genomeDir = REFERENCE_DIR + "genome/",
        sjdb = 100
    conda: "envs/star.yaml"
    threads: 24
    shell:
        """
        mkdir -p {params.genomeDir}; \
        gunzip -k {input.fasta}; \
        gunzip -k {input.gtf}; \
        FASTA="./reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa"; \
        GTF="./reference/Homo_sapiens.GRCh38.105.gtf"; \
        STAR --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeDir {params.genomeDir} \
        --genomeFastaFiles ${{FASTA}} \
        --sjdbGTFfile ${{GTF}} \
        --sjdbOverhang {params.sjdb} \
        """

rule star_map:
    input:
        fq1 = RESULT_DIR + "trimmed/" + "{sample}_R1_trimmed.fq.gz",
        fq2 = RESULT_DIR + "trimmed/" + "{sample}_R2_trimmed.fq.gz"
    output:
        RESULT_DIR + "star/{sample}_Aligned.sortedByCoord.out.bam",
        RESULT_DIR + "star/{sample}_Log.final.out"
    params:
        genomeDir = REFERENCE_DIR + "genome",
        prefix = RESULT_DIR + "star/{sample}_"
    conda: "envs/star.yaml"
    threads: 10
    shell:
        """
        STAR --genomeDir {params.genomeDir} \
        --runThreadN {threads} \
        --readFilesIn {input.fq1} {input.fq2} \
        --readFilesCommand zcat \
        --outFileNamePrefix {params.prefix} \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts
        """
        
rule generate_mapping_summary:
    input:
        expand(RESULT_DIR + "star/{sample}_Log.final.out", sample = SAMPLES)
    output:
        RESULT_DIR + "mapping_summary.csv"
    message:
        "Concatenating STAR mapping report and generating .csv mapping summary."
    params:
        directory_with_mapping_reports = RESULT_DIR + "star/",
        star_directory_name = RESULT_DIR + "star/"
    shell:
        "python scripts/generate_mapping_summary.py {params.directory_with_mapping_reports} {params.star_directory_name} {output}"
        
rule import_salmon_quant:
    input:
        expand(RESULT_DIR + "salmon_quant/" + "{sample}" + "/quant.sf", sample=SAMPLES)
    output:
        RESULT_DIR + "tximeta/" + "tximeta.rds"
    params:
        db = REFERENCE_DIR,
        metadata = config["samples_tbl"],
        salmon_dir = RESULT_DIR + "salmon_quant",
        salmon_index = REFERENCE_DIR + "salmon_index",
        fasta = config["refs"]["ensembl_transc_fasta"],
        gtf = config["refs"]["ensemblGTF"],
        jsonfile = REFERENCE_DIR + "annotations.json"
    conda: 
        "envs/tximeta.yaml"
    threads: 32
    shell:
        """
        Rscript --vanilla scripts/tximeta_import.R {params.db} {params.metadata} {params.salmon_dir} {output} {params.salmon_index} {params.fasta} {params.gtf} {params.jsonfile}
        """
        
rule gene_expression:
    input:
        RESULT_DIR + "tximeta/" + "tximeta.rds"
    output:
        dds = RESULT_DIR + "gene_expression/" + "dds.rds",
        dge = RESULT_DIR + "gene_expression/" + "dge.xlsx",
        counts = RESULT_DIR + "gene_expression/" + "counts.tsv"
    params:
        metadata = config["samples_tbl"],
        destDir = REFERENCE_DIR,
        designFormula = config["gexpression"]["design"]
    conda: "envs/deg.yaml"
    threads: 24
    resources: deseq=1
    shell:
        "Rscript --vanilla scripts/dge.R {input} {params.destDir} {params.metadata} {params.designFormula} {output.dds} {output.dge} {output.counts}"