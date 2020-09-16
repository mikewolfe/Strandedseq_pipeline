## HELPER FUNCTIONS inherited from parent SnakeFile:
# samples(pep)
# lookup_sample_metadata(sample, key, pep)
## TODO
# allow bwtools to deal with genomes with more than one chromosome or contig


def determine_resolution(config):
    try:
        resolution = config["coverage"]["resolution"]
    except KeyError:
        print(
        """
        Could not find specification for resolution in config file. I.e.

        coverage:
            resolution: N

        defaulting to 5 bp resolution
        """)
        resolution = 5
        
    return resolution

def determine_within_normalization(config):
    try:
        within = config["normalization"]["within"]
    except KeyError:
        print(
        """
        Could not find specification for normalization within a sample in config file. I.e.

        normalization:
            within: X

        defaulting to normalization by median signal in bins
        """)
        within = "median"
        
    return within

# SET CONSTANTS FOR COVERAGE CALCULATIONS
RES = determine_resolution(config)
WITHIN = determine_within_normalization(config)

def determine_files_for_log2ratio(config, pep, within):
    ending = "log2ratio"
    if "RobustZ" in config["normalization"]:
        RZ = config["normalization"]["RobustZ"]
        if RZ:
            ending += "RZ"
    samples = determine_extracted_samples(pep) 
    outfiles = ["results/deeptools_log2ratio/%s_%s_%s.bw"%(sample, within, ending) for sample in samples]
    return outfiles


# OVERALL RULES

rule clean_coverage_and_norm:
    shell:
        "rm -fr results/deeptools_coverage/ & "
        "rm -fr results/deeptools_log2ratio/ & "
        "rm -fr results/logs/deeptools_coverage/ & "
        "rm -fr results/logs/deeptools_log2ratio"

rule get_raw_coverage:
    input:
        expand("results/deeptools_coverage/{sample}_raw.bw", sample = samples(pep))

rule get_log2ratio_coverage:
    input:
        determine_files_for_log2ratio(config, pep, WITHIN)

def determine_genome_size(sample, config, pep):
    genome = lookup_sample_metadata(sample, "genome", pep) 
    try:
        genome_size = config["reference"][genome]["genome_size"]
    except KeyError:
        raise KeyError(
        """
        Cannot find genome_size for reference %s, make sure genome_size is
        specified in the config file. I.e.
        reference:
            %s:
                genome_size: NNNNN
        """%(genome, genome))
                
    return genome_size

rule deeptools_coverage_raw:
    input:
        inbam="results/bowtie2/{sample}_sorted.bam",
        ind="results/bowtie2/{sample}_sorted.bam.bai"
    output:
        "results/deeptools_coverage/{sample}_raw.bw"
    log:
        stdout="results/logs/deeptools_coverage/{sample}_raw.log",
        stderr="results/logs/deeptools_coverage/{sample}_raw.err"
    params:
        resolution = RES
    threads:
        5
    shell:
        "bamCoverage --bam {input.inbam} --outFileName {output} "
        "--outFileFormat 'bigwig' "
        "--numberOfProcessors {threads} --binSize {params.resolution} "
        "--samFlagInclude 66 --extendReads "
        "--exactScaling "
        "--minMappingQuality 10 > {log.stdout} 2> {log.stderr}"

rule deeptools_coverage:
    input:
        inbam="results/bowtie2/{sample}_sorted.bam",
        ind="results/bowtie2/{sample}_sorted.bam.bai"
    output:
        "results/deeptools_coverage/{sample}_{norm}.bw"
    log:
        stdout="results/logs/deeptools_coverage/{sample}_{norm}.log",
        stderr="results/logs/deeptools_coverage/{sample}_{norm}.err"
    params:
        resolution = RES,
        genome_size = lambda wildcards: determine_genome_size(wildcards.sample, config, pep)
    wildcard_constraints:
        norm="RPKM|CPM|BPM|RPGC"
    threads:
        5
    shell:
        "bamCoverage --bam {input.inbam} --outFileName {output} --outFileFormat 'bigwig' "
        "--numberOfProcessors {threads} --binSize {params.resolution} "
        "--samFlagInclude 66 --extendReads "
        "--normalizeUsing {wildcards.norm} "
        "--effectiveGenomeSize {params.genome_size} "
        "--exactScaling "
        "--minMappingQuality 10 > {log.stdout} 2> {log.stderr}"

# helper for finding correct input sample
def get_matching_input(sample, norm, pep):
   input_sample = lookup_sample_metadata(sample, "input_sample", pep)
   return "results/deeptools_coverage/%s_%s.bw"%(input_sample, norm)

rule deeptools_log2ratio:
    input:
        ext= "results/deeptools_coverage/{sample}_{norm}.bw",
        inp= lambda wildcards: get_matching_input(wildcards.sample,wildcards.norm, pep)
    output:
        "results/deeptools_log2ratio/{sample}_{norm}_log2ratio.bw"
    wildcard_constraints:
        norm="RPKM|CPM|BPM|RPGC|median"
    params: 
        resolution = RES
    log:
        stdout="results/logs/deeptools_log2ratio/{sample}_{norm}_log2ratio.log",
        stderr="results/logs/deeptools_log2ratio/{sample}_{norm}_log2ratio.err"
    threads:
        5
    shell:
        "bigwigCompare -b1 {input.ext} -b2 {input.inp} --outFileName {output} "
        "--operation 'log2' "
        "--binSize {params.resolution} "
        "--numberOfProcessors {threads} > {log.stdout} 2> {log.stderr}"

def get_input_bam(sample, pep):
    input_sample = lookup_sample_metadata(sample, "input_sample", pep)
    return "results/bowtie2/%s_sorted.bam"%(input_sample),

def get_input_bai(sample, pep):
    input_sample = lookup_sample_metadata(sample, "input_sample", pep)
    return "results/bowtie2/%s_sorted.bam.bai"%(input_sample),


rule deeptools_SES_log2ratio:
    input:
        ext= "results/bowtie2/{sample}_sorted.bam",
        ind_ext= "results/bowtie2/{sample}_sorted.bam.bai",
        inp = lambda wildcards: get_input_bam(wildcards.sample, pep),
        inp_ext = lambda wildcards: get_input_bai(wildcards.sample, pep)

    output:
        "results/deeptools_log2ratio/{sample}_SES_log2ratio.bw"
    log:
        stdout="results/logs/deeptools_log2ratio/{sample}_SES_log2ratio.log",
        stderr="results/logs/deeptools_log2ratio/{sample}_SES_log2ratio.err"
    threads:
        5
    params: 
        resolution = RES,
    shell:
        "bamCompare --bamfile1 {input.ext} --bamfile2 {input.inp} --outFileName {output} "
        "--outFileFormat 'bigwig' --scaleFactorsMethod 'SES' --operation 'log2' "
        "--extendReads --binSize {params.resolution} --samFlagInclude 66 "
        "--numberOfProcessors {threads} "
        "--exactScaling "
        "--minMappingQuality 10 > {log.stdout} 2> {log.stderr}"

rule bwtools_median:
    input:
        "results/deeptools_coverage/{sample}_raw.bw"
    output:
        "results/deeptools_coverage/{sample}_median.bw"
    params:
        resolution = RES,
        genome_size = lambda wildcards: determine_genome_size(wildcards.sample, config, pep),
        chrom_name = lambda wildcards: lookup_sample_metadata(wildcards.sample, "genome", pep)
    log:
        stdout="logs/bwtools/{sample}_median.log",
        stderr="logs/bwtools/{sample}_median.err"
    shell:
       "python3 "
       "workflow/scripts/bwtools.py "
       "{input} {output} --fr bigwig --to bigwig --chrm_name "
       "{params.chrom_name} --chrm_length {params.genome_size} "
       "--res {params.resolution} --operation Median_norm "
       "> {log.stdout} 2> {log.stderr}"

rule bwtools_RobustZ:
    input:
        "results/deeptools_log2ratio/{sample}_{norm}_log2ratio.bw"
    output:
        "results/deeptools_log2ratio/{sample}_{norm}_log2ratioRZ.bw"
    log:
        stdout="results/logs/bwtools/{sample}_{norm}_log2ratioRZ.log",
        stderr="results/logs/bwtools/{sample}_{norm}_log2ratioRZ.err"
    wildcard_constraints:
        norm="RPKM|CPM|BPM|RPGC|count|SES|median"
    params:
        resolution = RES,
        genome_size = lambda wildcards: determine_genome_size(wildcards.sample, config, pep),
        chrom_name = lambda wildcards: lookup_sample_metadata(wildcards.sample, "genome", pep)
    shell:
       "python3 "
       "workflow/scripts/bwtools.py "
       "{input} {output} --fr bigwig --to bigwig --chrm_name "
       "{params.chrom_name} "
       "--chrm_length {params.genome_size} " 
       "--res {params.resolution} "
       "--operation RobustZ > {log.stdout} 2> {log.stderr}"
