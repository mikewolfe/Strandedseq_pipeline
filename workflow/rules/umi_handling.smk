rule umi_extract_se:
    input:
        in1=lambda wildcards: match_fastq_to_sample(wildcards.sample, 'R0', pep, check_umi = False)
    output:
        out1=temp("results/umi_handling/umi_extract/{sample}_umi_extract_R0.fastq.gz"),
        out2="results/umi_handling/umi_extract/{sample}_umi_extract_R0_noumi.fastq.gz"
    params:
        umi_barcode_string = lambda wildcards: lookup_in_config_persample(config, pep, \
        ["preprocessing", "umi_extract_se", "umi_barcode_string"], wildcards.sample, \
        default = " "),
        umi_params = lambda wildcards: lookup_in_config_persample(config, pep, \
        ["preprocessing", "umi_extract_se", "umi_params"], wildcards.sample, \
        default = " ")
    log:
        stdout="results/umi_handling/logs/umi_extract_se/{sample}_umi_extract_se.log",
        stderr="results/umi_handling/logs/umi_extract_se/{sample}_umi_extract_se.err"
    conda:
        "../envs/umi_handling.yaml"
    shell:
        "umi_tools extract -I {input.in1:q} --bc-pattern={params.umi_barcode_string} "
        "-S {output.out1} --filtered-out={output.out2} {params.umi_params} > {log.stdout} 2> {log.stderr}"

rule umi_bam_sort:
    input: 
        inbam="results/alignment/bowtie2/unsorted/{sample}.bam"
    output:
        temp("results/umi_handling/bams/{sample}_sorted.bam")
    log:
        stderr="results/umi_handling/logs/bams/{sample}_sort.log"
    conda:
        "../envs/alignment.yaml"
    shell:
        "samtools sort {input} > {output} 2> {log.stderr}"

rule umi_bam_index:
    input:
        "results/umi_handling/bams/{sample}_sorted.bam"
    output:
        temp("results/umi_handling/bams/{sample}_sorted.bam.bai")
    log:
        stdout="results/umi_handling/logs/bams/{sample}_index.log",
        stderr="results/umi_handling/logs/bams/{sample}_index.log"
    conda:
        "../envs/alignment.yaml"
    shell:
        "samtools index {input} {output} > {log.stdout} 2> {log.stderr}"

rule umi_dedup:
    input:
        inbam="results/umi_handling/bams/{sample}_sorted.bam",
        indx="results/umi_handling/bams/{sample}_sorted.bam.bai"
    output:
        outbam=temp("results/umi_handling/umi_dedup/{sample}_dedup.bam")
    params:
        umi_params = lambda wildcards: lookup_in_config_persample(config, pep, \
        ["alignment", "umi_dedup", "umi_params"], wildcards.sample, \
        default = " "),
        stats="results/umi_handling/logs/umi_dedup/{sample}_umi_dedup_stats"
    log:
        stdout="results/umi_handling/logs/umi_dedup/{sample}_umi_dedup.log",
        stderr="results/umi_handling/logs/umi_dedup/{sample}_umi_dedup.err",
    conda:
        "../envs/umi_handling.yaml"
    shell:
        "umi_tools dedup -I {input.inbam} --output-stats={params.stats}  "
        "-S {output.outbam} > {log.stdout} 2> {log.stderr}"
