
def determine_NETseq_pause_output(config, pep):
    outfiles = []
    for model in lookup_in_config(config, ["NETseq", "pause_calling"], []):
        model_type = lookup_in_config(config, ["NETseq", "pause_calling", model, "model_type"], 
        err = "Need model type specified for pause_calling model %s in config file. I.e. \nNETseq:\n\tpause_calling:\n\t%s:\n\t\tmodel_type: 'Genomewide'"%(model,model))
        these_samples = filter_samples(pep, \
        lookup_in_config(config, ["NETseq", "pause_calling", model, "filter"], "not sample_name.isnull()"))
        if model_type == "Genomewide":
            for sample in these_samples:
                outfiles.append("results/NETseq/%s_pause/%s/%s_pause_calls.tsv.gz"%(model_type,model,sample))
        elif model_type == "Region":
            for sample in these_samples:
                outfiles.append("results/NETseq/%s_pause/%s/%s_pause_calls.tsv.gz"%(model_type,model,sample))
        else:
            raise ValueError("Model type %s not supported for NETseq pause calling. Need one of Genomewide or Region"%(model_type))
    return outfiles

def determine_logo_output(config, pep):
    out = []
    for val in determine_NETseq_pause_output(config, pep):
        out.append(val.replace("pause_calls.tsv.gz", "pause_logo.pdf"))
    return out
        

rule run_NETseq:
    input:
        determine_NETseq_pause_output(config, pep)

rule run_NETseq_logos:
    input:
        determine_logo_output(config, pep)

rule NETseq_genomewide_pause:
    input:
        "results/coverage_and_norm/deeptools_coverage/{sample}_plus_raw.bw",
        "results/coverage_and_norm/deeptools_coverage/{sample}_minus_raw.bw"
    output:
        "results/NETseq/Genomewide_pause/{model}/{sample}_pause_calls.tsv.gz"
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
        "results/NETseq/Region_pause/{model}/{sample}_pause_calls.tsv.gz"
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

def find_sample_fasta(sample, pep):
    genome = lookup_sample_metadata(sample, "genome", pep)
    return "results/alignment/combine_fasta/%s/%s.fa"%(genome, genome)

rule NETseq_pause_seqs:
    input:
        bed = "results/NETseq/{model_type}_pause/{model}/{sample}_pause_calls.tsv.gz",
        fasta = lambda wildcards: find_sample_fasta(wildcards.sample, pep)
    output:
        "results/NETseq/{model_type}_pause/{model}/{sample}_pause_seqs.txt.gz"
    log:
        stderr = "results/NETseq/logs/{model_type}_pause/{model}_{sample}_pause_seqs.err"
    params:
        upstream = 12,
        downstream = 5
    conda:
        "../envs/coverage_and_norm.yaml"
    threads:
        1
    shell:
        "zcat {input.bed} | python3 workflow/scripts/pull_seq_from_bed.py - {input.fasta} "
        "--upstream {params.upstream} --downstream {params.downstream} --circular "
        " 2> {log.stderr} | gzip >  {output}"

rule NETseq_pause_logo:
    input:
        "results/NETseq/{model_type}_pause/{model}/{sample}_pause_seqs.txt.gz"
    output:
        "results/NETseq/{model_type}_pause/{model}/{sample}_pause_logo.pdf"
    log:
        stderr = "results/NETseq/logs/{model_type}_pause/{model}_{sample}_pause_logo.err"
    threads:
        1
    conda:
        "../envs/R.yaml"
    shell:
        "Rscript workflow/scripts/NETseq_plot_logo.R {input} {output}"
