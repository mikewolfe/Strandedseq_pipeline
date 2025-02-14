
def determine_NETseq_pause_output(config, pep):
    outfiles = []
    for model in lookup_in_config(config, ["NETseq", "pause_calling"], []):
        model_type = lookup_in_config(config, ["NETseq", "pause_calling", model, "model_type"], 
        err = "Need model type specified for pause_calling model %s in config file. I.e. \nNETseq:\n\tpause_calling:\n\t%s:\n\t\tmodel_type: 'Genomewide'"%(model,model))
        get_seqs = lookup_in_config(config, ["NETseq", "pause_calling", model, "get_pause_seqs"], "true") == "true"
        these_samples = filter_samples(pep, \
        lookup_in_config(config, ["NETseq", "pause_calling", model, "filter"], "not sample_name.isnull()"))
        if model_type == "Genomewide":
            for sample in these_samples:
                outfiles.append("results/NETseq/%s_pause/%s/%s_pause_calls.bed.gz"%(model_type,model,sample))
                if get_seqs:
                    outfiles.append("results/NETseq/%s_pause/%s/%s_pause_seqs.txt.gz"%(model_type,model,sample))
        elif model_type == "Region":
            for sample in these_samples:
                outfiles.append("results/NETseq/%s_pause/%s/%s_pause_calls.bed.gz"%(model_type,model,sample))
                if get_seqs:
                    outfiles.append("results/NETseq/%s_pause/%s/%s_pause_seqs.txt.gz"%(model_type,model,sample))
        elif model_type == "Gini":
            for sample in these_samples:
                outfiles.append("results/NETseq/%s_pause/%s/%s_pause_calls.bed.gz"%(model_type,model,sample))
                if get_seqs:
                    outfiles.append("results/NETseq/%s_pause/%s/%s_pause_seqs.txt.gz"%(model_type,model,sample))
        else:
            raise ValueError("Model type %s not supported for NETseq pause calling. Need one of Genomewide or Region"%(model_type))
    return outfiles

def determine_logo_output(config, pep):
    outfiles = []
    for model in lookup_in_config(config, ["NETseq", "pause_calling"], []):
        if lookup_in_config(config, ["NETseq", "pause_calling", model, "plot_logos"], "true") != "true":
            continue
        model_type = lookup_in_config(config, ["NETseq", "pause_calling", model, "model_type"], 
        err = "Need model type specified for pause_calling model %s in config file. I.e. \nNETseq:\n\tpause_calling:\n\t%s:\n\t\tmodel_type: 'Genomewide'"%(model,model))
        these_samples = filter_samples(pep, \
        lookup_in_config(config, ["NETseq", "pause_calling", model, "filter"], "not sample_name.isnull()"))
        if model_type == "Genomewide":
            for sample in these_samples:
                outfiles.append("results/NETseq/%s_pause/%s/%s_pause_logo.pdf"%(model_type,model,sample))
        elif model_type == "Region":
            for sample in these_samples:
                outfiles.append("results/NETseq/%s_pause/%s/%s_pause_logo.pdf"%(model_type,model,sample))
        elif model_type == "Gini":
            for sample in these_samples:
                outfiles.append("results/NETseq/%s_pause/%s/%s_pause_logo.pdf"%(model_type,model,sample))
        else:
            raise ValueError("Model type %s not supported for NETseq pause calling. Need one of Genomewide or Region"%(model_type))
    return outfiles
        

rule run_NETseq:
    input:
        determine_NETseq_pause_output(config, pep)

rule run_NETseq_logos:
    input:
        determine_logo_output(config, pep)


rule NETseq_genomewide_gini:
    input:
        "results/coverage_and_norm/deeptools_coverage/{sample}_plus_raw.bw",
        "results/coverage_and_norm/deeptools_coverage/{sample}_minus_raw.bw"
    output:
        "results/NETseq/Gini_pause/{model}/{sample}_pause_calls.bed.gz"
    log:
        stderr = "results/NETseq/logs/Gini_pause/{sample}_gw_gini_{model}.err"
    params:
       NETseq_pause_params = lambda wildcards: lookup_in_config(config,
       ["NETseq", "pause_calling", wildcards.model, "NETseq_pause_params"],
       "--method Gini_null_ratio --n_boot 2000"),
       wsize = lambda wildcards: lookup_in_config(config,\
       ["NETseq", "pause_calling", wildcards.model, "wsize"],
       "100")
    conda:
        "../envs/coverage_and_norm.yaml"
    threads:
        10
    shell:
        "python3 "
        "workflow/scripts/NETseq_pause_calling.py Genomewide_gini {input} {params.wsize} "
        "{params.NETseq_pause_params} --p {threads} "
        "2> {log.stderr} | gzip > {output} "

rule NETseq_genomewide_pause:
    input:
        "results/coverage_and_norm/deeptools_coverage/{sample}_plus_raw.bw",
        "results/coverage_and_norm/deeptools_coverage/{sample}_minus_raw.bw"
    output:
        "results/NETseq/Genomewide_pause/{model}/{sample}_pause_calls.bed.gz"
    log:
        stderr = "results/NETseq/logs/Genomewide_pause/{sample}_gw_pause_{model}.err"
    params:
       NETseq_pause_params = lambda wildcards: lookup_in_config(config,
       ["NETseq", "pause_calling", wildcards.model, "NETseq_pause_params"],
       "--method Poisson --circular True"),
       wsize = lambda wildcards: lookup_in_config(config,\
       ["NETseq", "pause_calling", wildcards.model, "wsize"],
       "100")
    conda:
        "../envs/coverage_and_norm.yaml"
    threads:
        10
    shell:
        "python3 "
        "workflow/scripts/NETseq_pause_calling.py Genomewide {input} {params.wsize} "
        "{params.NETseq_pause_params} --p {threads} "
        "2> {log.stderr} | gzip > {output} "

rule NETseq_region_pause_calling:
    input:
        "results/coverage_and_norm/deeptools_coverage/{sample}_plus_raw.bw",
        "results/coverage_and_norm/deeptools_coverage/{sample}_minus_raw.bw",
        lambda wildcards: lookup_in_config(config,["NETseq", "pause_calling", wildcards.model, "regions"])

    output:
        "results/NETseq/Region_pause/{model}/{sample}_pause_calls.bed.gz"
    log:
        stderr = "results/NETseq/logs/Region_pause/{sample}_region_pause_{model}.err"
    params:
       NETseq_pause_params = lambda wildcards: lookup_in_config(config,\
       ["NETseq", "pause_calling", wildcards.model, "NETseq_pause_params"],
       "--method Poisson --circular True")
    conda:
        "../envs/coverage_and_norm.yaml"
    threads:
        10
    shell:
        "python3 "
        "workflow/scripts/NETseq_pause_calling.py Region {input} "
        "{params.NETseq_pause_params} --p {threads} "
        "2> {log.stderr} | gzip > {output} "

def find_sample_fasta(sample, pep, ending = ".fa"):
    genome = lookup_sample_metadata(sample, "genome", pep)
    return "results/alignment/combine_fasta/%s/%s%s"%(genome, genome, ending)

rule NETseq_pause_seqs:
    input:
        bed = "results/NETseq/{model_type}_pause/{model}/{sample}_pause_calls.bed.gz",
        fasta = lambda wildcards: find_sample_fasta(wildcards.sample, pep)
    output:
        "results/NETseq/{model_type}_pause/{model}/{sample}_pause_seqs.txt.gz"
    log:
        stderr = "results/NETseq/logs/{model_type}_pause/{model}_{sample}_pause_seqs.err"
    params:
       NETseq_logo_params = lambda wildcards: lookup_in_config(config,\
       ["NETseq", "pause_calling", wildcards.model, "NETseq_logo_params"],
       "--upstream 12 --downstream 5")
    conda:
        "../envs/coverage_and_norm.yaml"
    threads:
        1
    shell:
        "zcat {input.bed} | python3 workflow/scripts/pull_seq_from_bed.py - {input.fasta} "
        "{params.NETseq_logo_params} "
        " 2> {log.stderr} | gzip >  {output}"

rule generate_bg_model:
    input:
        "results/alignment/combine_fasta/{genome}/{genome}.fa"
    output:
        "results/alignment/combine_fasta/{genome}/{genome}_mm.txt"
    log:
        stderr = "results/alignment/logs/fasta-get-markov/{genome}.log"
    conda:
        "../envs/motif_calling.yaml"
    threads:
        1
    shell:
        "fasta-get-markov {input} > {output} 2> {log.stderr}"


rule NETseq_pause_logo:
    input:
        seqs="results/NETseq/{model_type}_pause/{model}/{sample}_pause_seqs.txt.gz",
        bg_model = lambda wildcards: find_sample_fasta(wildcards.sample, pep, ending = "_mm.txt")
    output:
        "results/NETseq/{model_type}_pause/{model}/{sample}_pause_logo.pdf",
        "results/NETseq/{model_type}_pause/{model}/{sample}_pause_logo_bg_corrected.pdf",
    params:
       NETseq_logo_display_params = lambda wildcards: lookup_in_config(config,\
       ["NETseq", "pause_calling", wildcards.model, "NETseq_logo_display_params"],
       "13 5")

    log:
        stdout = "results/NETseq/logs/{model_type}_pause/{model}_{sample}_pause_logo.log",
        stderr = "results/NETseq/logs/{model_type}_pause/{model}_{sample}_pause_logo.err"
    threads:
        1
    conda:
        "../envs/R.yaml"
    shell:
        "Rscript workflow/scripts/NETseq_plot_logo.R {input.seqs} results/NETseq/{wildcards.model_type}_pause/{wildcards.model}/{wildcards.sample}_pause_logo "
        "{input.bg_model} {params.NETseq_logo_display_params} > {log.stdout}  2> {log.stderr}"
