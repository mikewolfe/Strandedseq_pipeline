## HELPER FUNCTIONS
def samples(pep):
    """
    Get all of the unique sample names
    """
    return pep.sample_table["sample_name"]

def lookup_sample_metadata(sample, key, pep):
    """
    Get sample metadata by key
    """
    return pep.sample_table.at[sample, key]

# General read QC
rule read_qc:
    input:
        expand("results/fastqc_raw/{sample}_{pair}_fastqc.html", sample=samples(pep), pair=["R1", "R2"]),
        expand("results/fastqc_processed/{sample}_trim_paired_{pair}_fastqc.html", sample=samples(pep), pair=["R1", "R2"])

def match_fastq_to_sample(sample, pair, pep):
    out = lookup_sample_metadata(sample, "file_path", pep)
    if pair == "R1":
        out + lookup_sample_metadata(sample, "filenameR1", pep)
    elif pair == "R2":
        out + lookup_sample_metadata(sample, "filenameR2", pep)
    else:
        raise ValueError("Pair must be R1 or R2 not %s"%pair)
    return out

rule fastqc_raw:
    message: "Running fastqc on {wildcards.basename}"
    input:
       lambda wildcards: match_fastq_to_sample(wildcards.sample, wildcards.pair, pep) 
    output:
        "results/fastqc_raw/{sample}_{pair}_fastqc.html"
    threads: 1
    resources:
        mem_mb=5000
    log:
        stdout="results/logs/fastqc_raw/{sample}_{pair}_raw.log",
        stderr="results/logs/fastqc_raw/{sample}_{pair}_raw.err"
    shell:
        "fastqc {input} -o results/fastqc_raw > {log.stdout} 2> {log.stderr}"

rule fastqc_processed:
    message: "Running fastqc on {wildcards.sample} {wildcards.pair}"
    input:
        "results/trimmomatic/{sample}_trim_paired_{pair}.fastq.gz"
    output:
        "results/fastqc_processed/{sample}_trim_paired_{pair}_fastqc.html"
    threads: 1
    resources:
        mem_mb=5000
    log:
        stdout="results/logs/fastqc_processed/{sample}_trim_paired_{pair}.log",
        stderr="results/logs/fastqc_processed/{sample}_trim_paired_{pair}.err"
    shell:
        "fastqc {input} -o results/fastqc_processed > {log.stdout} 2> {log.stderr}"

# ChIP QC will be done by groups thus we need to write a few helpers

def determine_grouping_category(config):
    if "group_by" in config["quality_control"]:
        grouping_category = config["quality_control"]["group_by"]
    else:
        grouping_category = "genome"
    return grouping_category

def get_sample_by_group(group, pep):
    samp_table = pep.sample_table
    grouping_category = determine_grouping_category(config)
    samples = samp_table.loc[samp_table[grouping_category] == group, "sample_name"]
    return samples.tolist()

def qc_groups(config, pep):
    column = determine_grouping_category(config)
    return pep.sample_table[column].unique().tolist()

# OVERALL RULE FOR ChIP QC


rule ChIP_QC:
    input:
        expand("results/deeptools_QC/fingerprint/{group}_fingerprint.png", group = qc_groups(config, pep)),
        #expand("deeptools_QC/{samps}_GCBias_freqs.txt", samps=samples(pep)),
        expand("results/deeptools_QC/corHeatmap/{group}_corHeatmap.png", group = qc_groups(config, pep)),
        expand("results/deeptools_QC/corScatter/{group}_corScatterplot.png", group = qc_groups(config, pep)),
        expand("results/deeptools_QC/fragment_sizes/{group}_bamPEFragmentSize.png", group = qc_groups(config, pep)),
        expand("results/deeptools_QC/plotCoverage_rmdups/{group}_plotCoverage_rmdups.png", group = qc_groups(config, pep)),
        expand("results/deeptools_QC/plotCoverage/{group}_plotCoverage.png", group = qc_groups(config, pep)),
        expand("results/deeptools_QC/PCA/{group}_plotPCA.png", group = qc_groups(config, pep))


rule deeptools_QC_fingerprint:
    input:
        lambda wildcards: ["results/bowtie2/%s_sorted.bam"%samp for samp in get_sample_by_group(wildcards.group, pep)],
        lambda wildcards: ["results/bowtie2/%s_sorted.bam.bai"%samp for samp in get_sample_by_group(wildcards.group, pep)]
    output:
        outplot="results/deeptools_QC/fingerprint/{group}_fingerprint.png",
        rawcounts="results/deeptools_QC/fingerprint/{group}_fingerprint_counts.txt",
        qualmetrics="results/deeptools_QC/fingerprint/{group}_fingerprint_qual_metrics.txt"
    params:
        labels = lambda wildcards: " ".join([samp for samp in get_sample_by_group(wildcards.group, pep)]),
        bams = lambda wildcards: " ".join(["results/bowtie2/%s_sorted.bam"%samp for samp in get_sample_by_group(wildcards.group, pep)])
    threads:
        10
    log:
        stdout="results/logs/deeptools_QC/fingerprint/{group}_fingerprint.log",
        stderr="results/logs/deeptools_QC/fingerprint/{group}_fingerprint.err"

    shell:
        "plotFingerprint -b {params.bams} --plotFile {output.outplot} "
        "--outRawCounts {output.rawcounts} --outQualityMetrics {output.qualmetrics} "
        "--labels {params.labels} "
        "--binSize 500 "
        "--numberOfSamples 9200 --numberOfProcessors {threads} --extendReads "
        "--minMappingQuality 10 "
        "--samFlagInclude 66 > {log.stdout} 2> {log.stderr}"

rule deeptools_QC_multiBamSummary:
    input: 
        lambda wildcards: ["results/bowtie2/%s_sorted.bam"%samp for samp in get_sample_by_group(wildcards.group, pep)],
        lambda wildcards: ["results/bowtie2/%s_sorted.bam.bai"%samp for samp in get_sample_by_group(wildcards.group, pep)]
    output:
        "results/deeptools_QC/{group}_coverage_matrix.npz"
    params:
        labels = lambda wildcards: " ".join([samp for samp in get_sample_by_group(wildcards.group, pep)]),
        bams = lambda wildcards: " ".join(["results/bowtie2/%s_sorted.bam"%samp for samp in get_sample_by_group(wildcards.group, pep)])
    threads:
        10
    log:
        stdout="results/logs/deeptools_QC/{group}_multiBAMSummary.log",
        stderr="results/logs/deeptools_QC/{group}_multiBAMSummary.err"

    shell:
        "multiBamSummary bins -b {params.bams} -o {output} "
        "--labels {params.labels} "
        "--binSize 1000 "
        "--numberOfProcessors {threads} --extendReads "
        "--minMappingQuality 10 "
        "--samFlagInclude 66 > {log.stdout} 2> {log.stderr}"

rule deeptools_QC_corHeatmap:
    input:
        "results/deeptools_QC/{group}_coverage_matrix.npz"
    output:
        plot="results/deeptools_QC/corHeatmap/{group}_corHeatmap.png",
        cormat="results/deeptools_QC/corHeatmap/{group}_corvalues.txt"
    log:
        stdout="results/logs/deeptools_QC/corHeatmap/{group}_corHeatmap.log",
        stderr="results/logs/deeptools_QC/corHetamap/{group}_corHeatmap.err"
    shell:
        "plotCorrelation --corData {input} --corMethod 'spearman' "
        "--whatToPlot 'heatmap' --plotFile {output.plot} "
        "--colorMap 'viridis' "
        "--outFileCorMatrix {output.cormat} > {log.stdout} 2> {log.stderr}"

rule deeptools_QC_plotPCA:
    input:
        "results/deeptools_QC/{group}_coverage_matrix.npz"
    output:
        plot="results/deeptools_QC/PCA/{group}_plotPCA.png",
        pcaout="results/deeptools_QC/PCA/{group}_PCA_data.txt"
    log:
        stdout="results/logs/deeptools_QC/PCA/{group}_plotPCA.log",
        stderr="results/logs/deeptools_QC/PCA/{group}_plotPCA.err"
    shell:
        "plotPCA --corData {input} "
        "--plotFile {output.plot} "
        "--outFileNameData {output.pcaout} > {log.stdout} 2> {log.stderr}"

rule deeptools_QC_corScatterplot:
    input:
        "results/deeptools_QC/{group}_coverage_matrix.npz"
    output:
        plot="results/deeptools_QC/corScatter/{group}_corScatterplot.png"
    log:
        stdout="results/logs/deeptools_QC/corScatter/{group}_corScatterplot.log",
        stderr="results/logs/deeptools_QC/corScatter/{group}_corScatterplot.err"
    shell:
        "plotCorrelation --corData {input} --corMethod 'spearman' "
        "--whatToPlot 'scatterplot' --plotFile {output.plot} "
        "> {log.stdout} 2> {log.stderr}"

rule deeptools_QC_plotCoverage:
    input: 
        lambda wildcards: ["results/bowtie2/%s_sorted.bam"%samp for samp in get_sample_by_group(wildcards.group, pep)],
        lambda wildcards: ["results/bowtie2/%s_sorted.bam.bai"%samp for samp in get_sample_by_group(wildcards.group, pep)]
    output:
        plot="results/deeptools_QC/plotCoverage/{group}_plotCoverage.png",
        counts="results/deeptools_QC/plotCoverage/{group}_plotCoverage_count.txt"
    params:
        labels = lambda wildcards: " ".join([samp for samp in get_sample_by_group(wildcards.group, pep)]),
        bams = lambda wildcards: " ".join(["results/bowtie2/%s_sorted.bam"%samp for samp in get_sample_by_group(wildcards.group, pep)])
    threads:
        10
    log:
        stdout="results/logs/deeptools_QC/plotCoverage/{group}_plotCoverage.log",
        stderr="results/logs/deeptools_QC/plotCoverage/{group}_plotCoverage.err"

    shell:
        "plotCoverage -b {params.bams} -o {output.plot} "
        "--labels {params.labels} "
        "--numberOfProcessors {threads} --extendReads "
        "--numberOfSamples 1000 "
        "--outRawCounts {output.counts} "
        "--minMappingQuality 10 "
        "--samFlagInclude 66 > {log.stdout} 2> {log.stderr}"

rule deeptools_QC_plotCoverage_rmdups:
    input: 
        lambda wildcards: ["results/bowtie2/%s_sorted.bam"%samp for samp in get_sample_by_group(wildcards.group, pep)],
        lambda wildcards: ["results/bowtie2/%s_sorted.bam.bai"%samp for samp in get_sample_by_group(wildcards.group, pep)]
    output:
        plot="results/deeptools_QC/plotCoverage_rmdups/{group}_plotCoverage_rmdups.png",
        counts="results/deeptools_QC/plotCoverage_rmdups/{group}_plotCoverage_count_rmdups.txt"
    params:
        labels = lambda wildcards: " ".join([samp for samp in get_sample_by_group(wildcards.group, pep)]),
        bams = lambda wildcards: " ".join(["results/bowtie2/%s_sorted.bam"%samp for samp in get_sample_by_group(wildcards.group, pep)])
    threads:
        10
    log:
        stdout="results/logs/deeptools_QC/plotCoverage_rmdups/{group}_plotCoverage_rmdups.log",
        stderr="results/logs/deeptools_QC/plotCoverage_rmdups/{group}_plotCoverage_rmdups.err"
    shell:
        "plotCoverage -b {params.bams} --plotFile {output.plot} "
        "--labels {params.labels} "
        "--numberOfProcessors {threads} --extendReads --ignoreDuplicates "
        "--numberOfSamples 1000 "
        "--outRawCounts {output.counts} "
        "--minMappingQuality 10 "
        "--samFlagInclude 66 > {log.stdout} 2> {log.stderr}"

rule deeptools_QC_bamPEFragmentSize:
    input:
        lambda wildcards: ["results/bowtie2/%s_sorted.bam"%samp for samp in get_sample_by_group(wildcards.group, pep)],
        lambda wildcards: ["results/bowtie2/%s_sorted.bam.bai"%samp for samp in get_sample_by_group(wildcards.group, pep)]
    output:
        plot="results/deeptools_QC/fragment_sizes/{group}_bamPEFragmentSize.png",
        fragment_hist="results/deeptools_QC/fragment_sizes/{group}_bamPEFragmentSize.txt",
        table="results/deeptools_QC/fragment_sizes/{group}_bamPEFragmentSize_table.txt"
    params:
        labels = lambda wildcards: " ".join([samp for samp in get_sample_by_group(wildcards.group, pep)]),
        bams = lambda wildcards: " ".join(["results/bowtie2/%s_sorted.bam"%samp for samp in get_sample_by_group(wildcards.group, pep)])
    threads:
        10
    log:
        stdout="results/logs/deeptools_QC/fragment_sizes/{group}_bamPEFragmentSize.log",
        stderr="results/logs/deeptools_QC/fragment_sizes/{group}_bamPEFragmentSize.err"
    shell:
        "bamPEFragmentSize -b {params.bams} --histogram {output.plot} "
        "--samplesLabel {params.labels} "
        "--numberOfProcessors {threads} "
        "--distanceBetweenBins 10000 "
        "--table {output.table} --outRawFragmentLengths {output.fragment_hist}  > {log.stdout} 2> {log.stderr}"

# WORK IN PROGRESS
rule deeptools_QC_computeGCBias:
    input: 
        "results/bowtie2/{sample}_sorted.bam"
    output:
        gcfile="results/deeptools_QC/GCbias/{sample}_GCBias_freqs.txt",
        gcplot="results/deeptools_QC/GCbias/{sample}_GCBias_plot.png"
    threads:
        5
    log:
        stdout="results/logs/deeptools_QC/{sample}_GCBias.log",
        stderr="results/logs/deeptools_QC/{sample}_GCBias.err"
    params:
        genome_size = lambda wildcards: determine_genome_size(sample, config, pep)
    shell:
        "computeGCBias -b {input} --effectiveGenomeSize {params.genome_size} "
        "-g /mnt/scratch/mbwolfe/genomes/U00096_3.2bit "
        "--sampleSize 50000 "
        "--numberOfProcessors {threads} "
        "--GCbiasFrequenciesFile {output.gcfile}  > {log.stdout} 2> {log.stderr}"


