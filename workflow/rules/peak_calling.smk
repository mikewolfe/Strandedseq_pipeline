## HELPER FUNCTIONS inherited from parent SnakeFile:
# samples(pep)
# lookup_sample_metadata(sample, key, pep)

# GLOBAL CONSTANTS inherited from parent Snakefile:
# RES - resolution of coverage
# WITHIN - normalization to perform within samples
## TODO
# allow CMARRT to deal with genomes with more than one chromosome or contig

# HELPER FUNCTIONS
def determine_final_normalization(config, pep):
    ending = "log2ratio"
    if "normalization" in config and "RobustZ" in config["normalization"]:
        RZ = config["normalization"]["RobustZ"]
        if RZ:
            ending += "RZ"
    return ending
        
rule clean_peak_calling:
    shell:
        "rm -fr results/peak_calling"

rule run_peak_calling:
    input:
        expand("results/peak_calling/cmarrt/{sample}_{within}_{ending}.narrowPeak",\
        sample = determine_extracted_samples(pep),\
        within = WITHIN,\
        ending = determine_final_normalization(config, pep)),
        expand("results/peak_calling/macs2/{sample}_peaks.xls",\
        sample = determine_extracted_samples(pep))

# CMARRT helper functions
def determine_cmarrt_grouping(config):
    if "peak_calling" in config and \
       "cmarrt" in config["peak_calling"] and \
       "group_by" in config["peak_calling"]["cmarrt"]:
        grouping_category = config["peak_calling"]["cmarrt"]["group_by"]
    else:
        grouping_category = "genome"
    return grouping_category

def determine_cmarrt_params(sample, config, pep):
    samp_table = pep.sample_table
    grouping_category = determine_cmarrt_grouping(config)
    group = lookup_sample_metadata(sample, grouping_category, pep)
    if "peak_calling" in config and \
        "cmarrt" in config["peak_calling"] and \
        group in config["peak_calling"]["cmarrt"]:
        vals = config["peak_calling"]["cmarrt"][group]
    else:
        vals = {}
    if "wi" not in vals:
        vals["wi"] = 25
        logger.warning("CMARRT parameter wi not found for %s %s. Using default of %s"%(grouping_category, group, vals["wi"]))
    if "consolidate" not in vals:
        vals["consolidate"] = 10
        logger.warning("CMARRT parameter consolidate not found for %s %s. Using default of %s"%(grouping_category, group, vals["consolidate"]))    
    return vals

rule cmarrt_call_peaks:
    input:
        "results/coverage_and_norm/deeptools_log2ratio/{sample}_{norm}_{logratio}.npy"
    output:
        "results/peak_calling/cmarrt/{sample}_{norm}_{logratio}.narrowPeak"
    log:
        stdout="results/peak_calling/logs/cmarrt/{sample}_{norm}_{logratio}_cmarrt.log",
        stderr="results/peak_calling/logs/cmarrt/{sample}_{norm}_{logratio}_cmarrt.err"
    wildcard_constraints:
        norm="RPKM|CPM|BPM|RPGC|count|SES|median"
    params:
        peak_vals = lambda wildcards: determine_cmarrt_params(wildcards.sample, config, pep),
        outpre = lambda wildcards: "results/peak_calling/cmarrt/%s_%s_%s"%(wildcards.sample, wildcards.norm, wildcards.logratio),
        resolution = RES,
        genome_size = lambda wildcards: determine_genome_size(wildcards.sample, config, pep),
        chrom_name = lambda wildcards: lookup_sample_metadata(wildcards.sample, "genome", pep)
    conda:
        "../envs/peak_calling.yaml"
    shell:
        "python3 "
        "workflow/CMARRT_python/run_cmarrt.py "
        " {input} {params.peak_vals[wi]} "
        "{params.chrom_name} -o {params.outpre} "
        "--resolution {params.resolution} --np_start 0 "
        "--np_end {params.genome_size} --input_numpy "
        "--consolidate {params.peak_vals[consolidate]} --plots > {log.stdout} 2> {log.stderr}"


def get_macs2_matching_input(sample, pep):
   input_sample = lookup_sample_metadata(sample, "input_sample", pep)
   return "results/alignment/bowtie2/%s_sorted.bam"%(input_sample)

def determine_macs2_params(sample, config, pep):
    samp_table = pep.sample_table
    grouping_category = determine_cmarrt_grouping(config)
    group = lookup_sample_metadata(sample, grouping_category, pep)
    if "peak_calling" in config and \
        "macs2" in config["peak_calling"] and \
        group in config["peak_calling"]["macs2"]:
        vals = config["peak_calling"]["macs2"][group]
    else:
        vals = {}
    if "broad" not in vals:
        vals["broad"] = True
        logger.warning("Macs2 parameter broad not found for %s %s. Using default of %s"%(grouping_category, group, vals["broad"]))
    if vals["broad"]:
        vals["broad"] = "--broad"
    else:
        vals["broad"] = ""
    return vals

rule macs2_call_peaks:
    input:
        ext = "results/alignment/bowtie2/{sample}_sorted.bam",
        inp = lambda wildcards: get_macs2_matching_input(wildcards.sample, pep)
    output:
        "results/peak_calling/macs2/{sample}_peaks.xls"
    log:
        stdout="results/peak_calling/logs/macs2/{sample}_macs2.log",
        stderr="results/peak_calling/logs/macs2/{sample}_macs2.err"
    params:
        genome_size = lambda wildcards: determine_genome_size(wildcards.sample, config, pep),
        peak_vals = lambda wildcards: determine_macs2_params(wildcards.sample, config, pep),
        
    conda:
        "../envs/peak_calling.yaml"
    shell:
        "macs2 callpeak -t {input.ext} -c {input.inp} -n {wildcards.sample} "
        "--outdir results/peak_calling/macs2/ "
        "-f BAMPE  {params.peak_vals[broad]} "
        "-g {params.genome_size} > {log.stdout} 2> {log.stderr}"

