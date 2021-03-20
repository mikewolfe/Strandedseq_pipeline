
# OVERALL RULES

rule clean_postprocessing:
    shell:
        "rm -fr results/postprocessing/"

rule run_postprocessing:
    input:       
        expand("results/postprocessing/deeptools_byregion/{model}_deeptools_byregion.tab",\
        model = config["postprocessing"]["deeptools_byregion"].keys()),
        expand("results/postprocessing/deeptools_scaledregion/{model}_deeptools_scaledregion.tab",\
        model = config["postprocessing"]["deeptools_scaledregion"].keys()),
        expand("results/postprocessing/deeptools_referencepoint/{model}_deeptools_referencepoint.tab",\
        model = config["postprocessing"]["deeptools_referencepoint"].keys()),



def pull_bws_for_deeptools_models(toolname, modelname, config, pep):
    these_samples = filter_samples(pep, config["postprocessing"][toolname][modelname]["filter"])
    files = [config["postprocessing"][toolname][modelname]["filesignature"]%(sample) for sample in these_samples]
    return files

def pull_labels_for_deeptools_models(toolname, modelname, config, pep):
    these_samples = filter_samples(pep, config["postprocessing"][toolname][modelname]["filter"])
    return " ".join(these_samples)

def pull_param_for_postprocessing_model(toolname, modelname, param, config):
    return config["postprocessing"][toolname][modelname][param]

rule deeptools_byregion:
    input:
        inbws= lambda wildcards: pull_bws_for_deeptools_models("deeptools_byregion",wildcards.model,config, pep),
        inbed= lambda wildcards: pull_param_for_postprocessing_model("deeptools_byregion", wildcards.model, "regions", config)
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
        inbed= lambda wildcards: pull_param_for_postprocessing_model("deeptools_scaledregion", wildcards.model, "regions", config)
    output:
        outbinary="results/postprocessing/deeptools_scaledregion/{model}_deeptools_scaledregion.npz",
        outtext="results/postprocessing/deeptools_scaledregion/{model}_deeptools_scaledregion.tab"
    log:
        stdout="results/postprocessing/logs/deeptools_scaledregion/{model}.log",
        stderr="results/postprocessing/logs/deeptools_scaledregion/{model}.err"
    params:
        labels = lambda wildcards: pull_labels_for_deeptools_models("deeptools_scaledregion", wildcards.model, config, pep),
        upstream = lambda wildcards: pull_param_for_postprocessing_model("deeptools_scaledregion", wildcards.model, "upstream", config),
        downstream = lambda wildcards: pull_param_for_postprocessing_model("deeptools_scaledregion", wildcards.model, "downstream", config),
        scaleto = lambda wildcards: pull_param_for_postprocessing_model("deeptools_scaledregion", wildcards.model, "scaleto", config),
        binsize = lambda wildcards: pull_param_for_postprocessing_model("deeptools_scaledregion", wildcards.model, "binsize", config),
        binoperation = lambda wildcards: pull_param_for_postprocessing_model("deeptools_scaledregion", wildcards.model, "binoperation", config)
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
        inbed= lambda wildcards: pull_param_for_postprocessing_model("deeptools_referencepoint", wildcards.model, "regions", config)
    output:
        outbinary="results/postprocessing/deeptools_referencepoint/{model}_deeptools_referencepoint.npz",
        outtext="results/postprocessing/deeptools_referencepoint/{model}_deeptools_referencepoint.tab"
    log:
        stdout="results/postprocessing/logs/deeptools_referencepoint/{model}.log",
        stderr="results/postprocessing/logs/deeptools_referencepoint/{model}.err"
    params:
        labels = lambda wildcards: pull_labels_for_deeptools_models("deeptools_referencepoint", wildcards.model, config, pep),
        upstream = lambda wildcards: pull_param_for_postprocessing_model("deeptools_referencepoint", wildcards.model, "upstream", config),
        downstream = lambda wildcards: pull_param_for_postprocessing_model("deeptools_referencepoint", wildcards.model, "downstream", config),
        referencepoint = lambda wildcards: pull_param_for_postprocessing_model("deeptools_referencepoint", wildcards.model, "referencepoint", config),
        binsize = lambda wildcards: pull_param_for_postprocessing_model("deeptools_referencepoint", wildcards.model, "binsize", config),
        binoperation = lambda wildcards: pull_param_for_postprocessing_model("deeptools_referencepoint", wildcards.model, "binoperation", config)
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
