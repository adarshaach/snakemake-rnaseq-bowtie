# ----------------------------------------------------- #
# EXAMPLE WORKFLOW                                      #
# ----------------------------------------------------- #

# module to fetch genome from NCBI or Ensemble
# -----------------------------------------------------
rule get_genome:
    output:
        path=directory("results/get_genome"),
        fasta="results/get_genome/genome.fasta",
        gff="results/get_genome/genome.gff",
    conda:
        "../envs/get_genome.yml"
    message:
        """--- Parsing genome GFF and FASTA files."""
    params:
        database=config["get_genome"]["database"],
        assembly=config["get_genome"]["assembly"],
        fasta=config["get_genome"]["fasta"],
        gff=config["get_genome"]["gff"],
    log:
        path="results/get_genome/log/get_genome.log",
    script:
        "../scripts/get_genome.py"


# module simulate reads
# -----------------------------------------------------
rule simulate_reads:
    input:
        fasta=rules.get_genome.output.fasta,
    output:
        fastq1="results/simulate_reads/{sample}.bwa.read1.fastq.gz",
        fastq2="results/simulate_reads/{sample}.bwa.read2.fastq.gz",
    conda:
        "../envs/simulate_reads.yml"
    message:
        """--- Simulating reads with DWGSIM."""
    params:
        read_length=config["simulate_reads"]["read_length"],
        read_number=config["simulate_reads"]["read_number"],
        random_freq=config["simulate_reads"]["random_freq"],
    log:
        path="results/simulate_reads/{sample}.log",
    shell:
        "prefix=`echo {output.fastq1} | cut -f 1 -d .`;"
        "dwgsim -1 {params.read_length} -2 {params.read_length} "
        "-N {params.read_number} -o 1 -y {params.random_freq} {input.fasta} ${{prefix}} &> {log.path}"


rule simulate_umis:
    input:
        fastq1="results/simulate_reads/{sample}.bwa.read1.fastq.gz",
        fastq2="results/simulate_reads/{sample}.bwa.read2.fastq.gz",
    output:
        fastq1_umi="results/simulate_umis/{sample}.bwa.read1_umi.fastq.gz",
        fastq2_umi="results/simulate_umis/{sample}.bwa.read2_umi.fastq.gz",
    conda:
        "../envs/simulate_reads.yml"
    message:
        """--- Simulating reads with umis."""
    log:
        path="results/simulate_umis/{sample}_umi.log",
    shell:
        "zcat {input.fastq1} | awk '{{if(NR%4==1) print $1\":UMI:\"substr(\"ACGT\", int(rand()*4)+1, 1)substr(\"ACGT\", int(rand()*4)+1, 1)substr(\"ACGT\", int(rand()*4)+1, 1)substr(\"ACGT\", int(rand()*4)+1, 1); else print $0}}' | gzip > {output.fastq1_umi}; "
        "zcat {input.fastq2} | awk '{{if(NR%4==1) print $1\":UMI:\"substr(\"ACGT\", int(rand()*4)+1, 1)substr(\"ACGT\", int(rand()*4)+1, 1)substr(\"ACGT\", int(rand()*4)+1, 1)substr(\"ACGT\", int(rand()*4)+1, 1); else print $0}}' | gzip > {output.fastq2_umi};"

# module to make QC report
# -----------------------------------------------------
rule fastqc:
    input:
        fastq="results/simulate_umis/{sample}.bwa.{read}_umi.fastq.gz",
        
    output:
        report="results/fastqc/{sample}.bwa.{read}_umi_fastqc.html",
        
    conda:
        "../envs/fastqc.yml"
    message:
        """--- Checking fastq files with FastQC."""
    log:
        path="results/fastqc/{sample}.bwa.{read}.log",
    shell:
        "prefix=`echo {output.report} | cut -f 1-2 -d /`;"
        "fastqc --nogroup --extract --quiet --threads {threads} -o ${{prefix}} {input.fastq} > {log}"


# module to trim adapters from reads
# -----------------------------------------------------
rule cutadapt:
    input:
        fastq1="results/simulate_umis/{sample}.bwa.read1_umi.fastq.gz",
        fastq2="results/simulate_umis/{sample}.bwa.read2_umi.fastq.gz",
        
    output:
        fastq1="results/cutadapt/{sample}.bwa.read1.fastq.gz",
        fastq2="results/cutadapt/{sample}.bwa.read2.fastq.gz",
        
    conda:
        "../envs/cutadapt.yml"
    message:
        """--- Trim adapters from reads."""
    params:
        threep_adapter=config["cutadapt"]["threep_adapter"],
        fivep_adapter=config["cutadapt"]["fivep_adapter"],
        default=config["cutadapt"]["default"],
    log:
        stdout="results/cutadapt/{sample}.bwa.log",
        stderr="results/cutadapt/{sample}.bwa.stderr",
    threads: int(workflow.cores * 0.25)
    shell:
        "cutadapt {params.threep_adapter} {params.fivep_adapter} {params.default} --cores {threads} "
        "-o {output.fastq1} -p {output.fastq2} {input.fastq1} {input.fastq1} > {log.stdout} 2> {log.stderr}"
        

rule umi_extraction:
    input:
          fastq="results/cutadapt/{sample}.bwa.{read}.fastq.gz",
        
    output:
          fastq="results/umi_extraction/{sample}.bwa.{read}.fastq.gz",
        
    conda:
        "../envs/umitools.yml"
    message:
        """--- Extracting UMIs."""
    params:
        method=config["umi_extraction"]["method"],
        pattern=config["umi_extraction"]["pattern"],
    log:
        path="results/umi_extraction/log/{sample}.bwa.{read}.log",
    shell:
        "umi_tools extract --extract-method={params.method} --bc-pattern='{params.pattern}' --stdin {input.fastq} --stdout {output.fastq} > {log.path}"


rule bowtie_index:
    input:
        rules.get_genome.output.fasta,
    output:
        index=multiext(
            "results/bowtie_index/genome",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    conda:
        "../envs/mapping.yml"
    params:
        prefix=lambda wildcards, input: os.path.splitext(input[0])[0].replace(
            "results/get_genome", "results/bowtie_index"
        ),
    log:
        stdout="results/bowtie_index/bowtie2.log",
        stderr="results/bowtie_index/bowtie2.stderr",
    shell:
        "bowtie2-build {input} {params.prefix} > {log.stdout} 2> {log.stderr}"


rule mapping:
    input:
        r1="results/cutadapt/{sample}.bwa.read1.fastq.gz",
        r2="results/cutadapt/{sample}.bwa.read2.fastq.gz",
        index="results/bowtie_index/genome.1.bt2",
    output:
        bam="results/mapping/{sample}_mapping.bam",
        sorted_bam="results/mapping/sorted/{sample}_sorted.bam",
        bai="results/mapping/sorted/{sample}_sorted.bam.bai",
    log:
        path="results/mapping/{sample}_mapping.log",
    conda:
        "../envs/mapping.yml"
    threads: workflow.cores
    shell:
        "idx=`echo {input.index} | cut -f 1 -d .`;"
        "bowtie2 -x ${{idx}} -1 {input.r1} -2 {input.r2} -S {output.bam} --threads {threads} 2> {log.path};"
        "samtools sort -o {output.sorted_bam} {output.bam};"
        "samtools index {output.sorted_bam} {output.bai};"


rule umi_dedup:
    input:
        bam="results/mapping/sorted/{sample}_sorted.bam",
        bai="results/mapping/sorted/{sample}_sorted.bam.bai",
    output:
        dedup_bam="results/deduplicated/{sample}_dedup.bam",
        
    conda:
        "../envs/umitools.yml"
    params:
        default=config["umi_dedup"],
    log:
        path="results/deduplicated/log/{sample}.log",
        stderr="results/deduplicated/log/{sample}.stderr",
        stats="results/deduplicated/log/{sample}_umi_stats.txt",
       
    shell:
        "umi_tools dedup {params.default} -I {input.bam} -S {output.dedup_bam}"

# module to run multiQC on input + processed files
# -----------------------------------------------------
rule multiqc:
    input:
        expand(
            "results/simulate_reads/{sample}.bwa.{read}.fastq.gz",
            sample=samples.index,
            read=["read1", "read2"],
        ),
        expand(
            "results/simulate_umis/{sample}.bwa.{read}_umi.fastq.gz",
            sample=samples.index,
            read=["read1", "read2"],
        ),
        "results/bowtie_index/genome.1.bt2",
        expand(
            "results/cutadapt/{sample}.bwa.{read}.fastq.gz",
            sample=samples.index,
            read=["read1", "read2"],
        ),
        expand(
            "results/umi_extraction/{sample}.bwa.{read}.fastq.gz",
            sample=samples.index,
            read=["read1", "read2"],
        ),
        expand(
            "results/fastqc/{sample}.bwa.{read}_umi_fastqc.html",
            sample=samples.index,
            read=["read1", "read2"],
        ),
        expand(
            "results/mapping/{sample}_mapping.bam",
            sample=samples.index,
        ),
        expand(
            "results/mapping/sorted/{sample}_sorted.bam",
            sample=samples.index,
        ),
        expand(
            "results/deduplicated/{sample}_dedup.bam",
            sample=samples.index,
        ),
    output:
        report="results/multiqc/multiqc_report.html",
    conda:
        "../envs/multiqc.yml"
    message:
        """--- Generating MultiQC report for seq data."""
    params:
        config=config["multiqc"]["config"],
    log:
        path="results/multiqc/log/multiqc.log",
    shell:
        "outdir=`echo {output.report} | cut -f 1-2 -d /`; "
        "multiqc -c {params.config} --force --verbose --dirs --outdir ${{outdir}} results &> {log.path}"
