## HELPER FUNCTIONS inherited from parent SnakeFile:
# samples(pep)
# lookup_sample_metadata(sample, key, pep)

def get_fnas_for_bt2_index(config, reference):
    return config["reference"][reference]["fastas"]

def get_bt2_index(sample, pep):
    reference = lookup_sample_metadata(sample, "genome", pep)
    return "results/bowtie2_indices/%s/%s"%(reference,reference)

def get_bt2_index_file(sample, pep):
    reference = lookup_sample_metadata(sample, "genome", pep)
    return "results/bowtie2_indices/%s/%s.1.bt2"%(reference,reference)
    

rule pull_genbank:
    message: "Download genbank for {wildcards.accession}"
    output:
        "resources/genbanks/{accession}.gbk"
    log:
        stdout="results/logs/pull_genbank/{accession}.log",
        stderr="results/logs/pull_genbank/{accession}.err"
    threads: 1
    shell:
        "ncbi-acc-download {wildcards.accession} --out {output} > {log.stdout} "
        "2> {log.stderr}"

rule process_genbank:
    message: "Processing genbank for {wildcards.genome}"
    input:
        "resources/genbanks/{genome}.gbk"
    output:
        "results/process_genbank/{genome}/{genome}.{outfmt}"
    log:
        stdout="results/logs/process_genbank/{genome}_{outfmt}.log",
        stderr="results/logs/process_genbank/{genome}_{outfmt}.err"
    threads: 1
    shell:
         "python3 workflow/scripts/parse_genbank.py {input} "
         "--outfmt {wildcards.outfmt} "
         "--chrm '{wildcards.genome}'  "
         " > {output} 2> {log.stderr}"

rule bowtie2_index:
    input:
        lambda wildcards: get_fnas_for_bt2_index(config, wildcards.reference)
    output:
        "results/bowtie2_indices/{reference}/{reference}.1.bt2"
    params:
        fastas = lambda wildcards, input: ",".join(input)
    threads:
        5
    log:
        stdout="results/logs/bowtie2_index/{reference}.log",
        stderr="results/logs/bowtie2_index/{reference}.err" 
    shell:
        "bowtie2-build --threads {threads} "
        "{params.fastas} "
        "results/bowtie2_indices/{wildcards.reference}/{wildcards.reference} "
        "> {log.stdout} 2> {log.stderr}"

rule bowtie2_map:
    input:
        in1="results/trimmomatic/{sample}_trim_paired_R1.fastq.gz",
        in2="results/trimmomatic/{sample}_trim_paired_R2.fastq.gz",
        bt2_index= lambda wildcards: get_bt2_index_file(wildcards.sample,pep)
    output:
        temp("results/bowtie2/{sample}.bam")
    log:
        stderr="results/logs/bowtie2/{sample}_bt2.log" 
    params:
        bt2_index= lambda wildcards: get_bt2_index(wildcards.sample,pep)    
    threads: 
        5
    shell:
        "bowtie2 -x {params.bt2_index} -p {threads} "
        "-1 {input.in1} -2 {input.in2} --phred33  "
        "--end-to-end --very-sensitive 2> {log.stderr} "
        "| samtools view -b > {output}"

rule bam_sort:
    input:
        "results/bowtie2/{sample}.bam"
    output:
        "results/bowtie2/{sample}_sorted.bam"
    log:
        stderr="results/logs/bowtie2/{sample}_bt2_sort.log"
    shell:
        "samtools sort {input} > {output} 2> {log.stderr}"

rule bam_index:
    input:
        "results/bowtie2/{sample}_sorted.bam"
    output:
        "results/bowtie2/{sample}_sorted.bam.bai"
    log:
        stdout="results/logs/bowtie2/{sample}_bt2_index.log",
        stderr="results/logs/bowtie2/{sample}_bt2_index.log"
    shell:
        "samtools index {input} {output} > {log.stdout} 2> {log.stderr}"
