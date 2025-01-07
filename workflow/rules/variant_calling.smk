def which_samples_to_run(config, pep):
    these_samples = filter_samples(pep,
    lookup_in_config(config, ["variant_calling", "filter"], "not sample_name.isnull()"))
    return ["results/variant_calling/breseq/renamed_output/%s.gd"%sample for sample in these_samples]

rule run_variant_calling:
    input:
        "results/variant_calling/breseq/compiled_output.html"

rule clean_variant_calling:
    shell:
        "rm -fr results/variant_calling"
        

    
def get_references_per_sample(sample, pep, files = "genbanks_only"):
    '''
    Determine which files to pass as references to breseq.
    We want to preferentially use genbanks if possible but if
    no genbank exists than we can add the others as fastas.

    We also don't want to double add a fasta and a genbank for the
    same reference. This is brittle code to check for that case and
    only add fastas if they don't have a genbank file.

    For now we will run things assuming you only want to run on
    the genbanks by default
    '''
    ref_name = lookup_in_config_persample(config, pep, ["variant_calling", "reference"],
    sample, default = lookup_sample_metadata(sample, "genome", pep))
    this_genome = lookup_in_config(config, ["reference", ref_name])
    all_genbanks = []
    all_fastas = []
    basenames = set()
    if "genbanks" in this_genome:
        for key, value in this_genome["genbanks"].items():
            all_genbanks.append(value)
            basenames.add(key)
    if "fastas" in this_genome:
        for this_fasta in this_genome["fastas"]:
            all_fastas.append(this_fasta)

    # Logic to handle which reference files you want to include. For
    # now we will leave the default to only be genbanks
    if files == "genbanks_only":
        all_files = all_genbanks
    elif files == "fastas_only":
        all_files = all_fastas
    else:
    # deal with the tricky case of mixing fastas and genbanks
        all_files = all_genbanks
        for this_fasta in all_fastas:
            this_basename = os.path.splitext(os.path.basename(this_fasta))[0]
            if this_basename not in basenames:
                all_files.append(this_fasta)
                basenames.add(this_basename)
    return all_files

def get_all_refs(pep, files = "genbanks_only"): 
    these_samples = filter_samples(pep,
    lookup_in_config(config, ["variant_calling", "filter"], "not sample_name.isnull()"))
    out_files = set()
    for sample in these_samples:
        [out_files.add(this_file) for this_file in get_references_per_sample(sample, pep, files)]
    return list(out_files)

def format_references(all_files):
    out_str = ""
    for this_file in all_files:
        out_str += "-r %s "%(this_file)
    return out_str


def processed_fastqs(sample, pep):
    se = determine_single_end(sample, pep)
    out = []
    if se:
        out.append("results/preprocessing/trimmomatic/%s_trim_R0.fastq.gz"%(sample))
    else:
        out.append("results/preprocessing/trimmomatic/%s_trim_paired_R1.fastq.gz"%(sample))
        out.append("results/preprocessing/trimmomatic/%s_trim_paired_R2.fastq.gz"%(sample))
    return out

rule breseq:
    message: "Running breseq on {wildcards.sample}"
    input:
        processed_fastqs = lambda wildcards: processed_fastqs(wildcards.sample, pep),
        reference_files = lambda wildcards: get_references_per_sample(wildcards.sample, pep)
    output:
        "results/variant_calling/breseq/{sample}/output/summary.html",
        "results/variant_calling/breseq/{sample}/output/output.vcf",
        "results/variant_calling/breseq/{sample}/output/output.gd"
    threads: 10
    params:
        reference_file_string = lambda wildcards: format_references(get_references_per_sample(wildcards.sample, pep))
    log:
        stdout="results/variant_calling/logs/breseq/{sample}.log",
        stderr="results/variant_calling/logs/breseq/{sample}.err"
    conda:
        "../envs/variant_calling.yaml"
    shell:
        "breseq {params.reference_file_string} {input.processed_fastqs} "
        "-n {wildcards.sample} "
        "-o results/variant_calling/breseq/{wildcards.sample} "
        "-j {threads} > {log.stdout} 2> {log.stderr}"

rule rename_breseq_output:
    input:
        vcf="results/variant_calling/breseq/{sample}/output/output.vcf",
        gd="results/variant_calling/breseq/{sample}/output/output.gd"
    output:
        outvcf="results/variant_calling/breseq/renamed_output/{sample}.vcf",
        outgd="results/variant_calling/breseq/renamed_output/{sample}.gd"

    threads: 1
    shell:
        "cp {input.vcf} {output.outvcf} && cp {input.gd} {output.outgd}"

rule compile_breseq_output_html:
    input:
       which_samples_to_run(config, pep)
    output:
        "results/variant_calling/breseq/compiled_output.html"
    conda:
        "../envs/variant_calling.yaml"
    params:
        reference_file_string = format_references(get_all_refs(pep))
    log:
        stdout="results/variant_calling/logs/breseq/compiled_html.log",
        stderr="results/variant_calling/logs/breseq/compiled_html.err"
    shell:
        "gdtools ANNOTATE {params.reference_file_string} -f HTML "
        "-o {output} {input} > {log.stdout} 2> {log.stderr}"
        
