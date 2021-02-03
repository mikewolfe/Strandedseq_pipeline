## HELPER FUNCTIONS inherited from parent SnakeFile:
# samples(pep)
# lookup_sample_metadata(sample, key, pep)
# GLOBAL CONSTANTS inherited from parent Snakefile:
# RES - resolution of coverage
# WITHIN - normalization to perform within samples
# ENDING - log2ratio or log2ratioRZ


# OVERALL RULES

rule clean_coverage_and_norm:
    shell:
        "rm -fr results/coverage_and_norm/"

rule run_coverage_and_norm:
    input:       
        expand("results/coverage_and_norm/deeptools_log2ratio/{sample}_{within}_{ending}.bw",\
        sample = determine_extracted_samples(pep),\
        within = WITHIN,\
        ending = ENDING)

rule get_raw_coverage:
    input:
        expand("results/coverage_and_norm/deeptools_coverage/{sample}_raw.bw", sample = samples(pep))


def masked_regions_file_for_deeptools(config, sample, pep):
    genome = lookup_sample_metadata(sample, "genome", pep) 
    outfile = determine_masked_regions_file(config, genome)
    if outfile is not None:
        out = "--blackListFileName %s"%(outfile)
    else:
        out = ""
    return out

rule deeptools_coverage_raw:
    input:
        inbam="results/alignment/bowtie2/{sample}_sorted.bam",
        ind="results/alignment/bowtie2/{sample}_sorted.bam.bai"
    output:
        "results/coverage_and_norm/deeptools_coverage/{sample}_raw.bw"
    log:
        stdout="results/coverage_and_norm/logs/deeptools_coverage/{sample}_raw.log",
        stderr="results/coverage_and_norm/logs/deeptools_coverage/{sample}_raw.err"
    params:
        resolution = RES,
        masked_regions = lambda wildcards: masked_regions_file_for_deeptools(config, wildcards.sample, pep)
    threads:
        5
    conda:
        "../envs/coverage_and_norm.yaml"
    shell:
        "bamCoverage --bam {input.inbam} --outFileName {output} "
        "--outFileFormat 'bigwig' "
        "--numberOfProcessors {threads} --binSize {params.resolution} "
        "--samFlagInclude 66 --extendReads "
        "--exactScaling "
        "{params.masked_regions} "
        "--minMappingQuality 10 > {log.stdout} 2> {log.stderr}"

rule deeptools_coverage:
    input:
        inbam="results/alignment/bowtie2/{sample}_sorted.bam",
        ind="results/alignment/bowtie2/{sample}_sorted.bam.bai",
        genome_size= lambda wildcards: determine_effective_genome_size_file(wildcards.sample, config, pep)
    output:
        "results/coverage_and_norm/deeptools_coverage/{sample}_{norm}.bw"
    log:
        stdout="results/coverage_and_norm/logs/deeptools_coverage/{sample}_{norm}.log",
        stderr="results/coverage_and_norm/logs/deeptools_coverage/{sample}_{norm}.err"
    params:
        resolution = RES,
        genome_size = lambda wildcards: determine_effective_genome_size(wildcards.sample, config, pep),
        masked_regions = lambda wildcards: masked_regions_file_for_deeptools(config, wildcards.sample)
    wildcard_constraints:
        norm="RPKM|CPM|BPM|RPGC"
    threads:
        5
    conda:
        "../envs/coverage_and_norm.yaml"
    shell:
        "bamCoverage --bam {input.inbam} --outFileName {output} --outFileFormat 'bigwig' "
        "--numberOfProcessors {threads} --binSize {params.resolution} "
        "--samFlagInclude 66 --extendReads "
        "--normalizeUsing {wildcards.norm} "
        "--effectiveGenomeSize {params.genome_size} "
        "--exactScaling "
        "{params.masked_regions} "
        "--minMappingQuality 10 > {log.stdout} 2> {log.stderr}"

# helper for finding correct input sample
def get_log2ratio_matching_input(sample, norm, pep):
   input_sample = lookup_sample_metadata(sample, "input_sample", pep)
   return "results/coverage_and_norm/deeptools_coverage/%s_%s.bw"%(input_sample, norm)

rule deeptools_log2ratio:
    input:
        ext= "results/coverage_and_norm/deeptools_coverage/{sample}_{norm}.bw",
        inp= lambda wildcards: get_log2ratio_matching_input(wildcards.sample,wildcards.norm, pep)
    output:
        "results/coverage_and_norm/deeptools_log2ratio/{sample}_{norm}_log2ratio.bw"
    wildcard_constraints:
        norm="RPKM|CPM|BPM|RPGC|median"
    params: 
        resolution = RES
    log:
        stdout="results/coverage_and_norm/logs/deeptools_log2ratio/{sample}_{norm}_log2ratio.log",
        stderr="results/coverage_and_norm/logs/deeptools_log2ratio/{sample}_{norm}_log2ratio.err"
    threads:
        5
    conda:
        "../envs/coverage_and_norm.yaml"
    shell:
        "bigwigCompare -b1 {input.ext} -b2 {input.inp} --outFileName {output} "
        "--operation 'log2' "
        "--binSize {params.resolution} "
        "--numberOfProcessors {threads} > {log.stdout} 2> {log.stderr}"

def get_input_bam(sample, pep):
    input_sample = lookup_sample_metadata(sample, "input_sample", pep)
    return "results/alignment/bowtie2/%s_sorted.bam"%(input_sample),

def get_input_bai(sample, pep):
    input_sample = lookup_sample_metadata(sample, "input_sample", pep)
    return "results/alignment/bowtie2/%s_sorted.bam.bai"%(input_sample),


rule deeptools_SES_log2ratio:
    input:
        ext= "results/alignment/bowtie2/{sample}_sorted.bam",
        ind_ext= "results/alignment/bowtie2/{sample}_sorted.bam.bai",
        inp = lambda wildcards: get_input_bam(wildcards.sample, pep),
        inp_ext = lambda wildcards: get_input_bai(wildcards.sample, pep)

    output:
        "results/coverage_and_norm/deeptools_log2ratio/{sample}_SES_log2ratio.bw"
    log:
        stdout="results/coverage_and_norm/logs/deeptools_log2ratio/{sample}_SES_log2ratio.log",
        stderr="results/coverage_and_norm/logs/deeptools_log2ratio/{sample}_SES_log2ratio.err"
    threads:
        5
    params: 
        resolution = RES
    conda:
        "../envs/coverage_and_norm.yaml"
    shell:
        "bamCompare --bamfile1 {input.ext} --bamfile2 {input.inp} --outFileName {output} "
        "--outFileFormat 'bigwig' --scaleFactorsMethod 'SES' --operation 'log2' "
        "--extendReads --binSize {params.resolution} --samFlagInclude 66 "
        "--numberOfProcessors {threads} "
        "--exactScaling "
        "--minMappingQuality 10 > {log.stdout} 2> {log.stderr}"

rule bwtools_median:
    input:
        "results/coverage_and_norm/deeptools_coverage/{sample}_raw.bw"
    output:
        "results/coverage_and_norm/deeptools_coverage/{sample}_median.bw"
    params:
        resolution = RES
    log:
        stdout="results/coverage_and_norm/logs/bwtools/{sample}_median.log",
        stderr="results/coverage_and_norm/logs/bwtools/{sample}_median.err"
    conda:
        "../envs/coverage_and_norm.yaml"
    shell:
       "python3 "
       "workflow/scripts/bwtools.py "
       "{input} {output} "
       "--res {params.resolution} --operation Median_norm "
       "> {log.stdout} 2> {log.stderr}"

rule bwtools_RobustZ:
    input:
        "results/coverage_and_norm/deeptools_log2ratio/{sample}_{norm}_log2ratio.bw"
    output:
        "results/coverage_and_norm/deeptools_log2ratio/{sample}_{norm}_log2ratioRZ.bw"
    log:
        stdout="results/coverage_and_norm/logs/bwtools/{sample}_{norm}_log2ratioRZ.log",
        stderr="results/coverage_and_norm/logs/bwtools/{sample}_{norm}_log2ratioRZ.err"
    wildcard_constraints:
        norm="RPKM|CPM|BPM|RPGC|count|SES|median"
    params:
        resolution = RES
    conda:
        "../envs/coverage_and_norm.yaml"
    shell:
       "python3 "
       "workflow/scripts/bwtools.py "
       "{input} {output} "
       "--res {params.resolution} "
       "--operation RobustZ > {log.stdout} 2> {log.stderr}"
