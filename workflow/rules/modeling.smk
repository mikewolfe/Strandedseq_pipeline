rule clean_modeling:
    shell:
        "rm -rf results/modeling/"


def determine_modeling_output(config, pep):
    outstr = "results/modeling/%s/%s/%s.rds"
    outfiles = []
    for model in lookup_in_config(config, ['modeling'], {}).keys():
        program = lookup_in_config(config, ['modeling', model, 'program'], 'sleuth')
        outfiles.append(outstr%(model, program, model))
    return outfiles

rule run_modeling:
    input:
        determine_modeling_output(config, pep) 

def all_htseq_for_modeling(model, config, pep):
    these_samples = filter_samples(pep, \
    lookup_in_config(config, ["modeling", model, "filter"], "not sample_name.isnull()"))
    quant_path = lookup_in_config(config, ["modeling", model, "count_path"]) + "%s.tsv"
    files = [quant_path%(sample) for sample in these_samples]
    return files

def pull_sample_labels_for_modeling(modelname, config, pep):
    these_samples = filter_samples(pep, lookup_in_config(config, ["modeling", modelname, "filter"], "not sample_name.isnull()"))
    return " ".join(these_samples)

rule HTseq_readcount:
    input:
        inbam = "results/alignment/bowtie2/{sample}_sorted.bam",
        inbam_idx = "results/alignment/bowtie2/{sample}_sorted.bam.bai",
        ingff= lambda wildcards: lookup_in_config(config, ["feature_counting", wildcards.model, "regions"])
    output:
        outtext="results/modeling/feature_counting/{model}/HTseq_readcount/{sample}.tsv"
    log:
        stderr="results/modeling/logs/feature_counting/{model}/HTseq_readcount/{sample}_HTseq_readcount.err"
    params:
        HTseq_params = lambda wildcards: lookup_in_config(config,\
        ["feature_counting", wildcards.model, "HTseq_params"],\
        "-s no -a 0 -t ID")
    conda:
        "../envs/modeling.yaml"
    shell:
        "htseq-count -f bam -r pos "
        "{input.inbam} {input.ingff} "
        "{params.HTseq_params} " 
        "> {output.outtext} 2> {log.stderr} "

def create_metadata_table(model, config, pep, outfile):
    """
    create a metadata table from config file and kallisto paths
    """
    sample_table = pep.sample_table
    all_keys = list(sample_table.columns)
    all_keys.remove("sample_name")
    if "sample" in all_keys:
        all_keys.remove("sample")
    if "fileName" in all_keys:
        all_keys.remove("fileName")
    all_samples = filter_samples(pep, lookup_in_config(config, ["modeling", model, "filter"], "not sample_name.isnull()"))
    # write output to a file
    with open(outfile, mode = "w") as outf:
        # header needs sample and kallisto path
        header = "sample\tfileName\t" + "\t".join(all_keys) + "\n"
        outf.write(header)
        for sample in all_samples:
            # for any value that doesn't exist for that sample give an NA
            values = [str(lookup_sample_metadata(sample, key, pep)) for key in all_keys]
            line = "%s\t"%(sample) + "%s.tsv\t"%(sample)+ "\t".join(values) + "\n"
            outf.write(line)

rule DESeq2_metadata:
    message: "Creating metadata table for model {wildcards.model}"
    output:
        "results/modeling/{model}/DESeq2/{model}_metadata.tsv"
    threads: 1
    run:
        create_metadata_table(wildcards.model, config, pep, output[0])

rule DESeq2_diffexp:
    message: "Running differential expression for model {wildcards.model}"
    input: 
        metadata = "results/modeling/{model}/DESeq2/{model}_metadata.tsv",
        samples = lambda wildcards: all_htseq_for_modeling(wildcards.model, config, pep)
    output:
        model="results/modeling/{model}/DESeq2/{model}.rds",
        lrt="results/modeling/{model}/DESeq2/{model}_lrt.tsv",
        coefficients="results/modeling/{model}/DESeq2/{model}_coefficients.tsv",
        scaled_count="results/modeling/{model}/DESeq2/{model}_scaledcount.tsv"
    params:
        path_to_HTseq_files = lambda wildcards: lookup_in_config(config,\
        ['modeling', wildcards.model, 'count_path']),
        model_full = lambda wildcards: lookup_in_config(config,\
        ['modeling', wildcards.model, 'full']),
        model_reduced = lambda wildcards: lookup_in_config(config,\
        ['modeling', wildcards.model, 'reduced'],\
        "~ 1"),
        spike_regions = lambda wildcards: lookup_in_config(config,\
        ['modeling', wildcards.model, 'spike_regions'],\
        "NA")
    log:
        stdout="results/modeling/logs/{model}/DESeq2/{model}.log",
        stderr="results/modeling/logs/{model}/DESeq2/{model}.err"
    threads: 10 
    conda:
        "../envs/modeling.yaml"
    resources:
        mem_mb=5000
    shell:
        "Rscript workflow/scripts/DESeq2_diffexp.R {input.metadata} {params.path_to_HTseq_files} "
        "results/modeling/{wildcards.model}/DESeq2/{wildcards.model} "
        "'{params.model_full}' '{params.model_reduced}' {threads} {params.spike_regions}"
        " > {log.stdout} 2> {log.stderr}"


def get_all_kallisto_output_dirs(model, config):
    """
    Pull all the kallisto paths per sample
    """
    samples = filter_samples(pep, lookup_in_config(config, ["modeling", model, "filter"], "not sample_name.isnull()"))
    out_samples = {}
    for sample in samples:
        out_samples[sample] = "results/pseudoalignment/kallisto/%s/%s"%(model,sample)
    return out_samples

def get_all_kallisto_output_h5s(model, config):
    """
    Get all kallisto output h5s
    """
    sample_dict = get_all_kallisto_output_dirs(model, config)
    out_samples = [sample + "/abundance.h5" for sample in sample_dict.values()]
    return out_samples


def create_sleuth_metadata_table(model, config, pep, outfile):
    """
    create a metadata table from config file and kallisto paths
    """
    paths = get_all_kallisto_output_dirs(model, config)
    sample_table = pep.sample_table
    all_keys = list(sample_table.columns)
    all_keys.remove("sample_name")
    if "path" in all_keys:
        all_keys.remove("path")
    if "sample" in all_keys:
        all_keys.remove("path")
    all_samples = filter_samples(pep, lookup_in_config(config, ["modeling", model, "filter"], "not sample_name.isnull()"))
    # write output to a file
    with open(outfile, mode = "w") as outf:
        # header needs sample and kallisto path
        header = "sample\t" + "\t".join(all_keys) + "\tpath\n"
        outf.write(header)
        for sample in all_samples:
            # for any value that doesn't exist for that sample give an NA
            values = [str(lookup_sample_metadata(sample, key, pep)) for key in all_keys]
            line = "%s\t"%(sample) + "\t".join(values) + "\t" + paths[sample] + "\n"
            outf.write(line)


rule sleuth_metadata:
    message: "Creating metadata table for model {wildcards.model}"
    input:
        lambda wildcards: get_all_kallisto_output_h5s(wildcards.model, config)
    output:
        "results/modeling/{model}/sleuth/samples.tsv"
    threads: 1
    run:
        create_sleuth_metadata_table(wildcards.model, config, pep, output[0])

rule sleuth_diffexp:
    message: "Running differential expression for model {wildcards.model}"
    input: 
        "results/modeling/{model}/sleuth/samples.tsv"
    output:
        model="results/modeling/{model}/sleuth/{model}.rds",
        lrt="results/modeling/{model}/sleuth/{model}_lrt.tsv",
        coefficients="results/modeling/{model}/sleuth/{model}_coefficients.tsv",
        tpm="results/modeling/{model}/sleuth/{model}_normed_tpm.tsv"
    params:
        model_full = lambda wildcards: config['modeling'][wildcards.model]['full'],
        model_reduced = lambda wildcards: config['modeling'][wildcards.model]['reduced']
    log:
        stdout="results/modeling/logs/{model}/sleuth_diffexp/{model}.log",
        stderr="results/modeling/logs/{model}/sleuth_diffexp/{model}.err"
    threads: 10 
    conda:
        "../envs/modeling.yaml"
    resources:
        mem_mb=5000
    shell:
        "Rscript workflow/scripts/sleuth_diffexp.R {input} results/modeling/{wildcards.model}/sleuth/{wildcards.model} "
        "'{params.model_full}' '{params.model_reduced}' {threads} "
        "> {log.stdout} 2> {log.stderr}"
