rule clean_pseudoalignment:
   shell:
    "rm -rf results/pseudoalignment/"


rule bed_to_fasta:
    message: "Converting bed file into fasta file"
    input:
        to_convert= lambda wildcards: lookup_in_config(config, ["modeling", wildcards.model, "regions"]),
        genome=lambda wildcards: lookup_in_config(config, ["modeling", wildcards.model, "genome_fasta"])
    output:
        "results/pseudoalignment/transcriptome/{model}_transcriptome.fa"
    threads: 1
    conda:
        "../envs/alignment.yaml"
    log:
        stdout = "results/pseudoalignment/logs/bed_to_fasta/{model}.log",
        stderr = "results/pseudoalignment/logs/bed_to_fasta/{model}.err"
    shell:
        "python workflow/scripts/bed_to_fasta.py {input.to_convert} {input.genome} {output} > {log.stdout} 2> {log.stderr}"

rule kallisto_index:
    message: "Creating kallisto index for assembly {wildcards.model}"
    input:
        "results/pseudoalignment/transcriptome/{model}_transcriptome.fa"
    output:
        "results/pseudoalignment/kallisto/{model}/kallisto.idx"
    threads: 1
    conda:
        "../envs/pseudoalignment.yaml"
    log:
        stdout="results/pseudoalignment/logs/kallisto/{model}/kallisto_index.log",
        stderr="results/pseudoalignment/logs/kallisto/{model}/kallisto_index.err"
    shell:
        "kallisto index {input} -i {output} > {log.stdout} 2> {log.stderr}"

def configure_librarytype_for_kallisto(sample):
    library = lookup_sample_metadata(sample, 'library_type', pep)
    if library == "rf-stranded":
        out = "--rf-stranded"
    elif library == "fr-stranded":
        out = "--fr-stranded"
    elif "unstranded":
        out = ""
    else:
        raise ValueError("Only unstranded, rf-stranded, and fr-stranded are supported. Not %s"%(library))
    return out

def get_fastqs_for_kallisto(sample, pep):
    out = []
    if determine_single_end(sample, pep):
        out.append("results/preprocessing/trimmomatic/%s_trim_R0.fastq.gz"%sample)
    else:
        out.append("results/preprocessing/trimmomatic/%s_trim_paired_R1.fastq.gz"%sample)
        out.append("results/preprocessing/trimmomatic/%s_trim_paired_R2.fastq.gz"%sample)
    return out

def kallisto_single_end_params(config, pep, sample):
    if determine_single_end(sample, pep):

        out = lookup_in_config_persample(config, pep,["pseudoalignment", "kallisto", "single_end_params"],
        sample,
        default = "--single -l 200 -s 20")
    else:
        out = ""
    return out
    

rule kallisto_quant:
    message: "Quantifying {wildcards.sample} using transcriptome library type {params.library}, and params {params.kallisto_params}"
    input:
        index="results/pseudoalignment/kallisto/{model}/kallisto.idx",
        fastqs = lambda wildcards: get_fastqs_for_kallisto(wildcards.sample, pep)
    params:
        kallisto_params = lambda wildcards: lookup_in_config_persample(config, pep,
        ["pseudoalignment", "kallisto", "params"], wildcards.sample, default = "-b 100"),
        library = lambda wildcards: configure_librarytype_for_kallisto(wildcards.sample),
        singleend = lambda wildcards: kallisto_single_end_params(config, pep, wildcards.sample)
    threads: 10
    resources:
        mem_mb=10000
    output:
        out1="results/pseudoalignment/kallisto/{model}/{sample}/abundance.h5",
        out2="results/pseudoalignment/kallisto/{model}/{sample}/abundance.tsv",
    conda:
        "../envs/pseudoalignment.yaml"
    log:
        stdout="results/pseudoalignment/logs/kallisto/{model}/{sample}_align.log",
        stderr="results/pseudoalignment/logs/kallisto/{model}/{sample}_align.err"
    shell:
        "kallisto quant {input.fastqs} "
        "-i {input.index} -o results/pseudoalignment/kallisto/{wildcards.model}/{wildcards.sample} "
        "{params.kallisto_params} {params.library} {params.singleend} "
        "> {log.stdout} 2> {log.stderr} "
