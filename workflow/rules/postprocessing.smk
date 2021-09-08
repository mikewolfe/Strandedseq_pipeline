
# OVERALL RULES

rule clean_postprocessing:
    shell:
        "rm -fr results/postprocessing/"

def determine_postprocessing_files(config):
    outfiles = []
    if lookup_in_config(config, ["postprocessing", "deeptools_byregion"], ""):
        outfiles.extend(
        ["results/postprocessing/deeptools_byregion/%s_deeptools_byregion.tab"%model\
        for model in config["postprocessing"]["deeptools_byregion"]])

    if lookup_in_config(config, ["postprocessing", "deeptools_scaledregion"], ""):
        outfiles.extend(
        ["results/postprocessing/deeptools_byregion/%s_deeptools_scaledregion.tab"%model\
        for model in config["postprocessing"]["deeptools_scaledregion"]])


    if lookup_in_config(config, ["postprocessing", "deeptools_referencepoint"], ""):
        outfiles.extend(
        ["results/postprocessing/deeptools_byregion/%s_deeptools_referencepoint.tab"%model\
        for model in config["postprocessing"]["deeptools_referencepoint"]])


    if lookup_in_config(config, ["postprocessing", "bwtools_query"], ""):
        outfiles.extend(
        ["results/postprocessing/bwtools_query/%s_bwtools_query.tab"%model\
        for model in config["postprocessing"]["bwtools_query"]])
    return outfiles



rule run_postprocessing:
    input:       
        determine_postprocessing_files(config)



def pull_bws_for_deeptools_models(toolname, modelname, config, pep):
    these_samples = filter_samples(pep, \
    lookup_in_config(config, ["postprocessing", toolname, modelname, "filter"], "not input_sample.isnull()"))
    file_sig = lookup_in_config(config, ["postprocessing", toolname, modelname, "filesignature"],\
    "results/coverage_and_norm/bwtools_compare/%s_median_log2ratio.bw")
    files = [file_sig%(sample) for sample in these_samples]
    return files

def pull_labels_for_deeptools_models(toolname, modelname, config, pep):
    these_samples = filter_samples(pep, lookup_in_config(config, ["postprocessing", toolname, modelname, "filter"], "not input_sample.isnull()"))
    return " ".join(these_samples)

rule deeptools_byregion:
    input:
        inbws= lambda wildcards: pull_bws_for_deeptools_models("deeptools_byregion",wildcards.model,config, pep),
        inbed= lambda wildcards: lookup_in_config(config, ["postprocessing", "deeptools_byregion", wildcards.model, "regions"], None)
    output:
        outbinary="results/postprocessing/deeptools_byregion/{model}_deeptools_byregion.npz",
        outtext="results/postprocessing/deeptools_byregion/{model}_deeptools_byregion.tab"
    log:
        stdout="results/postprocessing/logs/deeptools_byregion/{model}.log",
        stderr="results/postprocessing/logs/deeptools_byregion/{model}.err"
    params:
        labels = lambda wildcards: pull_labels_for_deeptools_models("deeptools_byregion", wildcards.model, config, pep)
    threads:
        5
    conda:
        "../envs/coverage_and_norm.yaml"
    shell:
        "multiBigwigSummary BED-file "
        "--bwfiles {input.inbws} "
        "--outFileName {output.outbinary} "
        "--BED {input.inbed} "
        "--labels {params.labels} "
        "--outRawCounts {output.outtext} "
        "--numberOfProcessors {threads} > {log.stdout} 2> {log.stderr} "


rule deeptools_scaledregion:
    input:
        inbws= lambda wildcards: pull_bws_for_deeptools_models("deeptools_scaledregion",wildcards.model,config, pep),
        inbed= lambda wildcards: lookup_in_config(config, ["postprocessing", "deeptools_scaledregion", wildcards.model, "regions"], None)
    output:
        outbinary="results/postprocessing/deeptools_scaledregion/{model}_deeptools_scaledregion.npz",
        outtext="results/postprocessing/deeptools_scaledregion/{model}_deeptools_scaledregion.tab"
    log:
        stdout="results/postprocessing/logs/deeptools_scaledregion/{model}.log",
        stderr="results/postprocessing/logs/deeptools_scaledregion/{model}.err"
    params:
        labels = lambda wildcards: pull_labels_for_deeptools_models("deeptools_scaledregion", wildcards.model, config, pep),
        upstream = lambda wildcards: lookup_in_config(config, ["postprocessing", "deeptools_scaledregion", wildcards.model, "upstream"], 0),
        downstream = lambda wildcards: lookup_in_config(config, ["postprocessing", "deeptools_scaledregion", wildcards.model, "downstream"], 0),
        scaleto = lambda wildcards: lookup_in_config(config, ["postprocessing", "deeptools_scaledregion", wildcards.model, "scaleto"], 1000),
        binsize = lambda wildcards: lookup_in_config(config, ["postprocessing", "deeptools_scaledregion", wildcards.model, "binsize"], 5),
        binoperation = lambda wildcards: lookup_in_config(config, ["postprocessing", "deeptools_scaledregion", wildcards.model, "binoperation"], "mean")
    threads:
        5
    conda:
        "../envs/coverage_and_norm.yaml"
    shell:
        "computeMatrix scale-regions "
        "--scoreFileName {input.inbws} "
        "--outFileName {output.outbinary} "
        "--regionsFileName {input.inbed} "
        "--samplesLabel {params.labels} "
        "--outFileNameMatrix {output.outtext} "
        "--upstream {params.upstream} "
        "--downstream {params.downstream} "
        "--regionBodyLength {params.scaleto} "
        "--averageTypeBins {params.binoperation} "
        "--sortRegions keep "
        "--numberOfProcessors {threads} > {log.stdout} 2> {log.stderr} "

rule deeptools_referencepoint:
    input:
        inbws= lambda wildcards: pull_bws_for_deeptools_models("deeptools_referencepoint",wildcards.model,config, pep),
        inbed= lambda wildcards: lookup_in_config(config, ["postprocessing", "deeptools_referencepoint", wildcards.model, "regions"], None)
    output:
        outbinary="results/postprocessing/deeptools_referencepoint/{model}_deeptools_referencepoint.npz",
        outtext="results/postprocessing/deeptools_referencepoint/{model}_deeptools_referencepoint.tab"
    log:
        stdout="results/postprocessing/logs/deeptools_referencepoint/{model}.log",
        stderr="results/postprocessing/logs/deeptools_referencepoint/{model}.err"
    params:
        labels = lambda wildcards: pull_labels_for_deeptools_models("postprocessing", "deeptools_scaledregion", wildcards.model, config, pep),
        upstream = lambda wildcards: lookup_in_config(config, ["postprocessing", "deeptools_scaledregion", wildcards.model, "upstream"], 0),
        downstream = lambda wildcards: lookup_in_config(config, ["postprocessing", "deeptools_scaledregion", wildcards.model, "downstream"], 0),
        referencepoint = lambda wildcards: lookup_in_config(config, ["postprocessing", "deeptools_referencepoint", wildcards.model, "referencepoint"], "TSS"),
        binsize = lambda wildcards: lookup_in_config(config, ["postprocessing", "deeptools_scaledregion", wildcards.model, "binsize"], 5),
        binoperation = lambda wildcards: lookup_in_config(config, ["postprocessing", "deeptools_scaledregion", wildcards.model, "binoperation"], "mean")
    threads:
        5
    conda:
        "../envs/coverage_and_norm.yaml"
    shell:
        "computeMatrix reference-point "
        "--scoreFileName {input.inbws} "
        "--outFileName {output.outbinary} "
        "--regionsFileName {input.inbed} "
        "--samplesLabel {params.labels} "
        "--outFileNameMatrix {output.outtext} "
        "--upstream {params.upstream} "
        "--downstream {params.downstream} "
        "--referencePoint {params.referencepoint} "
        "--averageTypeBins {params.binoperation} "
        "--sortRegions keep "
        "--numberOfProcessors {threads} > {log.stdout} 2> {log.stderr} "


rule bwtools_query:
    input:
        inbws= lambda wildcards: pull_bws_for_deeptools_models("bwtools_query",wildcards.model,config, pep),
        inbed= lambda wildcards: lookup_in_config(config, ["postprocessing", "bwtools_query", wildcards.model, "regions"], None)
    output:
        outtext="results/postprocessing/bwtools_query/{model}_bwtools_query.tab"
    log:
        stdout="results/postprocessing/logs/bwtools_query/{model}.log",
        stderr="results/postprocessing/logs/bwtools_query/{model}.err"
    params:
        labels = lambda wildcards: pull_labels_for_deeptools_models("bwtools_query", wildcards.model, config, pep),
        upstream = lambda wildcards: lookup_in_config(config, ["postprocessing", "bwtools_query", wildcards.model, "upstream"], 0),
        downstream = lambda wildcards: lookup_in_config(config, ["postprocessing", "bwtools_query", wildcards.model, "downstream"], 0),
        res = lambda wildcards: lookup_in_config(config, ["postprocessing", "bwtools_query", wildcards.model, "res"], 5),
        summarize = lambda wildcards: lookup_in_config(config, ["postprocessing", "bwtools_query", wildcards.model, "summarize"], 'single'),
        summary_func = lambda wildcards: lookup_in_config(config, ["postprocessing", "bwtools_query", wildcards.model, "summary_func"], 'mean')
    threads:
        5
    conda:
        "../envs/coverage_and_norm.yaml"
    shell:
        "python3 workflow/scripts/bwtools.py query "
        "{output.outtext} "
        "{input.inbws} "
        "--res {params.res} "
        "--regions {input.inbed} "
        "--upstream {params.upstream} "
        "--downstream {params.downstream} "
        "--samp_names {params.labels} "
        "--summarize {params.summarize} "
        "--summary_func {params.summary_func} "
        "> {log.stdout} 2> {log.stderr} "
