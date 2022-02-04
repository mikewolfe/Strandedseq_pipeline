def get_bisulfite_samples():
    samples = filter_samples(pep, "seq_method == 'bisulfite'") 
    outfiles = ["results/bisulfite/bismark_extract/CpG_CTOT_%s_pe.deduplicated.txt.gz"%samp for samp in samples]
    return outfiles

rule run_bisulfite:
    input:
        get_bisulfite_samples()

rule bismark_genome_prep:
    message: "Generating bisulfite genome for genome {wildcards.genome}"
    input:
        "results/alignment/combine_fasta/{genome}/{genome}.fa"
    output:
        "results/alignment/combine_fasta/{genome}/Bisulfite_Genome/CT_conversion/BS_CT.1.bt2",
        "results/alignment/combine_fasta/{genome}/Bisulfite_Genome/GA_conversion/BS_GA.1.bt2"
    log:
        stdout="results/bisulfite/logs/bismark_genome_prep/{genome}/{genome}.log",
        stderr="results/bisulfite/logs/bismark_genome_prep/{genome}/{genome}.err"
    threads: 1
    params:
        path_to_genome = "results/alignment/combine_fasta/{genome}/"
    conda:
        "../envs/bisulfite.yaml"
    shell:
         "bismark_genome_preparation --verbose {params.path_to_genome} "
         "> {log.stdout} 2> {log.stderr}"

rule bismark_align:
    message: "Aligning {wildcards.sample} to bisulfite genome"
    input:
        in1="results/preprocessing/trimmomatic/{sample}_trim_paired_R1.fastq.gz",
        in2="results/preprocessing/trimmomatic/{sample}_trim_paired_R2.fastq.gz",
        genome= lambda wildcards: "results/alignment/combine_fasta/%s/Bisulfite_Genome/CT_conversion/BS_CT.1.bt2"%(\
                           lookup_sample_metadata(wildcards.sample, "genome", pep))
    output:
        "results/bisulfite/bismark_align/{sample}_pe.bam"
    log:
        stdout="results/bisulfite/logs/bismark_align/{sample}.log",
        stderr="results/bisulfite/logs/bismark_align/{sample}.err"
    threads: 5
    params:
        path_to_genome = lambda wildcards: "results/alignment/combine_fasta/%s/"%(lookup_sample_metadata(wildcards.sample, "genome", pep)),
        align_params = lambda wildcards: lookup_in_config_persample(config,\
        pep, ["bisulfite", "bismark_align", "bismark_align_param_string"], wildcards.sample,\
        "--pbat")
    conda:
        "../envs/bisulfite.yaml"
    shell:
         "bismark {params.path_to_genome} -1 {input.in1} -2 {input.in2} "
         "--output_dir results/bisulfite/bismark_align/ "
         "--basename {wildcards.sample} {params.align_params} "
         "> {log.stdout} 2> {log.stderr}"

rule bismark_dedup:
    input:
        "results/bisulfite/bismark_align/{sample}_pe.bam"
    output:
        "results/bisulfite/bismark_align/{sample}_pe.deduplicated.bam"
    log:
        stdout="results/bisulfite/logs/bismark_dedup/{sample}.log",
        stderr="results/bisulfite/logs/bismark_dedup/{sample}.err"
    threads: 1
    conda:
        "../envs/bisulfite.yaml"
    shell:
         "deduplicate_bismark  {input} --output_dir results/bisulfite/bismark_align/"
         "> {log.stdout} 2> {log.stderr}"

rule bismark_extract:
    input:
        "results/bisulfite/bismark_align/{sample}_pe.deduplicated.bam"
    output:
        "results/bisulfite/bismark_extract/CpG_CTOT_{sample}_pe.deduplicated.txt.gz",
        "results/bisulfite/bismark_extract/CpG_CTOB_{sample}_pe.deduplicated.txt.gz",
        "results/bisulfite/bismark_extract/Non_CpG_CTOT_{sample}_pe.deduplicated.txt.gz",
        "results/bisulfite/bismark_extract/Non_CpG_CTOB_{sample}_pe.deduplicated.txt.gz"
    log:
        stdout="results/bisulfite/logs/bismark_extract/{sample}.log",
        stderr="results/bisulfite/logs/bismark_extract/{sample}.err"
    threads: 1
    conda:
        "../envs/bisulfite.yaml"
    shell:
         "bismark_methylation_extractor  -o results/bisulfite/bismark_extract/ -p "
         "--no_header --no_overlap --merge_non_CpG {input} --gzip"
         "> {log.stdout} 2> {log.stderr}"
