from pathlib import Path

# control of the pipeline
configfile: "config/config.yaml"
# sample metadata and information
pepfile: "pep/config.yaml"

# include in several rules here
include: "workflow/rules/preprocessing.smk"
include: "workflow/rules/alignment.smk"
include: "workflow/rules/coverage_and_norm.smk"

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

def all_raw_fastqs(pep):
    """
    Get all of the raw fastq files
    """
    these_samples = samples(pep)
    read1s = [lookup_sample_metadata(sample, "file_path", pep) + 
    lookup_sample_metadata(sample, "filenameR1", pep) for sample in these_samples]
    read2s = [lookup_sample_metadata(sample, "file_path", pep) + 
    lookup_sample_metadata(sample, "filenameR2", pep) for sample in these_samples]
    read1s.extend(read2s)
    return read1s

def pull_processed_fastqs(samplename, config):
    """
    Figure out processed fastqname(s) from samplename
    """
    read1 = "results/trimmomatic/%s_trim_paired_R1.fastq.gz"%(samplename)
    read2 = "results/trimmomatic/%s_trim_paired_R2.fastq.gz"%(samplename)
    return([read1, read2])

## overall rules

rule clean:
    threads: 1
    shell:
        "rm -rf results/"
    
rule get_region_averages:
    input:
        expand("results/deeptools_regions/regionAverage_{gff}_{norm}.npz" , gff=config["gffs"], norm=["median_log2ratio_RZ"])

rule call_peaks:
    input:
        expand("results/cmarrt/{sampname}_{norm}.narrowPeak", sampname=config["named_ext"], norm=["median_log2ratio_RZ"])

rule summary_plots:
    input:
        expand("results/deeptools_summary/plotProfile_ref_{gff}_{norm}.png", gff=config["gffs"], norm=["median_log2ratio_RZ", "median_log2ratio"]),
        expand("results/deeptools_summary/plotProfile_scale_{gff}_{norm}.png", gff=config["gffs"], norm=["median_log2ratio_RZ", "median_log2ratio"]),
        expand("results/deeptools_summary/plotProfile_ref_{gff}_{norm}_heatmap.png", gff=config["gffs"], norm=["median_log2ratio_RZ", "median_log2ratio"])
#        expand("deeptools_summary/plotProfile_ref_{gff}_{norm}.png", gff=config["gffs"], norm=["SES_log2ratio_RZ", "BPM_log2ratio_RZ", "RPGC_log2ratio_RZ", "median_log2ratio_RZ", "count_log2ratio_RZ", "median_log2ratio"]),
#        expand("deeptools_summary/plotProfile_scale_{gff}_{norm}.png", gff=config["gffs"], norm=["SES_log2ratio_RZ", "BPM_log2ratio_RZ", "RPGC_log2ratio_RZ", "median_log2ratio_RZ", "count_log2ratio_RZ", "median_log2ratio"]),
        #expand("deeptools_summary/plotProfile_ref_{gff}_{norm}.png", gff=config["gffs"], norm=["SES_log2ratio", "BPM_log2ratio", "RPGC_log2ratio", "median_log2ratio", "count_log2ratio", "median_log2ratio"]),
        #expand("deeptools_summary/plotProfile_scale_{gff}_{norm}.png", gff=config["gffs"], norm=["SES_log2ratio", "BPM_log2ratio", "RPGC_log2ratio", "median_log2ratio", "count_log2ratio", "median_log2ratio"])

rule chIP_QC:
    input:
        expand("results/deeptools_QC/{inpname}_fingerprint.png", inpname=config["named_inp"]),
        #expand("deeptools_QC/{samps}_GCBias_freqs.txt", samps=config["samples"]),
        "results/deeptools_QC/all_corHeatmap.png",
        "results/deeptools_QC/all_corScatterplot.png",
        "results/deeptools_QC/all_bamPEFragmentSize.png",
        "results/deeptools_QC/all_plotCoverage_rmdups.png",
        "results/deeptools_QC/all_plotCoverage.png",
        "results/deeptools_QC/all_plotPCA.png"

rule get_logratios:
    input:
#        expand("deeptools_coverage/{sampname}_SES_log2ratio.bw", sampname=config["named_ext"]),
#        expand("deeptools_coverage/{sampname}_count_log2ratio.bw", sampname=config["named_ext"])
        expand("results/deeptools_coverage/{sampname}_median_log2ratio_RZ.bw", sampname=config["named_ext"])

rule get_raw_bedgraph:
    input:
        expand("results/raw_coverage/{sample}.bedgraph.gz", sample=config["samples"]),
        expand("results/deeptools_coverage/{sample}_raw.bedgraph", sample=config["samples"]),
        expand("results/deeptools_coverage/{sample}_RPGC.bedgraph", sample=config["samples"])

rule qc:
    input:
        expand("fastqc_raw/{sample}_{pair}_fastqc.html", sample=config["samples"], pair=["R1", "R2"]),
        expand("fastqc_processed/{sample}_trim_paired_{pair}_fastqc.html", sample=config["samples"], pair=["R1", "R2"])


rule fastqc_raw:
    message: "Running fastqc on {wildcards.basename}"
    input:
       lambda wildcards: match_fastq_to_basename(wildcards.basename) 
    output:
        "results/fastqc_raw/{basename}_fastqc.html"
    threads: 1
    resources:
        mem_mb=5000
    log:
        stdout="results/logs/fastqc_raw/{basename}_raw.log",
        stderr="results/logs/fastqc_raw/{basename}_raw.err"
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



rule bwtools_bw2npy:
    input:
        "results/deeptools_coverage/{sampname}_{norm}.bw"
    output:
        "results/deeptools_coverage/{sampname}_{norm}.npy"
    log:
        stdout="results/logs/bwtools/{sampname}_{norm}_bw2npy.log",
        stderr="results/logs/bwtools/{sampname}_{norm}_bw2npy.err"
    shell:
       "python3 "
       "workflow/scripts/bwtools.py "
       "{input} {output} --fr bigwig --to numpy --chrm_name "
       "'U00096.3' --chrm_length 4641652 --res 5 > {log.stdout} 2> {log.stderr}"

rule cmarrt_call_peaks:
    input:
        "deeptools_coverage/{sampname}_{norm}.npy"
    output:
        "cmarrt/{sampname}_{norm}.narrowPeak"
    log:
        stdout="logs/cmarrt/{sampname}_{norm}_cmarrt.log",
        stderr="logs/cmarrt/{sampname}_{norm}_cmarrt.err"
    params:
        wi=25,
        chrom_name="U00096.3",
        cons=10,
        outpre = lambda wildcards: "cmarrt/%s_%s"%(wildcards.sampname, wildcards.norm)
    shell:
        "python3 "
        "/mnt/Groups/LandickGrp/General/BobFile/USER_FILES/mbwolfe/tools/CMARRT_python/run_cmarrt.py "
        " {input} {params.wi} {params.chrom_name} -o {params.outpre} --resolution 5 --np_start 0 "
        "--np_end 4641652 --input_numpy --consolidate {params.cons} --plots > {log.stdout} 2> {log.stderr}"


rule deeptools_regionAverage:
    input:
       lambda wildcards: ["deeptools_coverage/%s_%s.bw"%(val, wildcards.norm) for val in config["named_ext"]]
    output:
        outfile="deeptools_regions/regionAverage_{gff}_{norm}.npz",
        outmat="deeptools_regions/regionAverage_{gff}_{norm}.tab"
    wildcard_constraints:
        norm="(RPKM|CPM|BPM|RPGC|median).*"
    log:
        stdout="logs/deeptools/{gff}_{norm}_regionAverage.log",
        stderr="logs/deeptools/{gff}_{norm}_regionAverage.err"
    params:        
        bwfiles = lambda wildcards: " ".join(
                    ["deeptools_coverage/%s_%s.bw"%(val,wildcards.norm) 
                    for val in config["named_ext"]]),
        gfffile = lambda wildcards: config["gffs"][wildcards.gff],
        names = " ".join([val for val in config["named_ext"]])
    threads:
        5
    shell:
        "multiBigwigSummary BED-file --bwfiles {params.bwfiles} "
        "--BED {params.gfffile} --outFileName {output.outfile} --outRawCounts {output.outmat} "
        "--labels {params.names} "
        "--numberOfProcessors {threads} > {log.stdout} 2> {log.stderr}"


rule deeptools_count_logratio:
    input:
        ext= lambda wildcards: "bowtie2/%s_sorted.bam"%(config["named_ext"][wildcards.sampname]),
        inp= lambda wildcards: "bowtie2/%s_sorted.bam"%(config["named_inp"][config["input_matching"][wildcards.sampname]]),
        ind_ext= lambda wildcards: "bowtie2/%s_sorted.bam"%(config["named_ext"][wildcards.sampname]),
        ind_inp= lambda wildcards: "bowtie2/%s_sorted.bam.bai"%(config["named_inp"][config["input_matching"][wildcards.sampname]])

    output:
        "deeptools_coverage/{sampname}_count_log2ratio.bw"
    log:
        stdout="logs/deeptools/{sampname}_count_log2ratio.log",
        stderr="logs/deeptools/{sampname}_count_log2ratio.err"
    threads:
        5
    shell:
        "bamCompare --bamfile1 {input.ext} --bamfile2 {input.inp} --outFileName {output} "
        "--outFileFormat 'bigwig' --scaleFactorsMethod 'readCount' --operation 'log2' "
        "--exactScaling --extendReads --binSize 5 --samFlagInclude 66 "
        "--numberOfProcessors {threads} "
        "--minMappingQuality 10 > {log.stdout} 2> {log.stderr}"

def get_input_for_fingerprint(wildcards):
    all_files = []
    for val in config["extracted_matching"][wildcards.inpname]:
        all_files.append("bowtie2/%s_sorted.bam"%config["named_ext"][val])
        all_files.append("bowtie2/%s_sorted.bam.bai"%config["named_ext"][val])
    all_files.append("bowtie2/%s_sorted.bam"%config["named_inp"][wildcards.inpname])
    all_files.append("bowtie2/%s_sorted.bam.bai"%config["named_inp"][wildcards.inpname])
    return all_files

rule deeptools_QC_fingerprint:
    input: get_input_for_fingerprint
    output:
        outplot="deeptools_QC/{inpname}_fingerprint.png",
        rawcounts="deeptools_QC/{inpname}_fingerprint_counts.txt",
        qualmetrics="deeptools_QC/{inpname}_fingerprint_qual_metrics.txt"
    params:
        contname=lambda wildcards: wildcards.inpname,
        contbams=lambda wildcards: "bowtie2/%s_sorted.bam"%config["named_inp"][wildcards.inpname],
        extname=lambda wildcards: " ".join([val for val in config["extracted_matching"][wildcards.inpname]]),
        extbams=lambda wildcards: " ".join(["bowtie2/%s_sorted.bam"%config["named_ext"][val] for val in config["extracted_matching"][wildcards.inpname]])

    threads:
        10
    log:
        stdout="logs/deeptools/{inpname}_fingerprint.log",
        stderr="logs/deeptools/{inpname}_fingerprint.err"

    shell:
        "plotFingerprint -b {params.contbams} {params.extbams} --plotFile {output.outplot} "
        "--outRawCounts {output.rawcounts} --outQualityMetrics {output.qualmetrics} "
        "--JSDsample {params.contbams} --labels {params.contname} {params.extname} "
        "--numberOfSamples 9200 --numberOfProcessors {threads} --extendReads "
        "--minMappingQuality 10 "
        "--samFlagInclude 66 > {log.stdout} 2> {log.stderr}"

rule deeptools_QC_multiBamSummary:
    input: expand("bowtie2/{sample}_sorted.bam", sample=[config["named_ext"][val] for val in config["named_ext"]]),
           expand("bowtie2/{sample}_sorted.bam.bai", sample=[config["named_ext"][val] for val in config["named_ext"]]),
           expand("bowtie2/{sample}_sorted.bam", sample=[config["named_inp"][val] for val in config["named_inp"]]),
           expand("bowtie2/{sample}_sorted.bam.bai", sample=[config["named_inp"][val] for val in config["named_inp"]]),
    output:
        "deeptools_QC/all_coverage_matrix.npz"
    params:
        extbams = " ".join(["bowtie2/%s_sorted.bam"%config["named_ext"][val] for val in config["named_ext"]]),
        inpbams = " ".join(["bowtie2/%s_sorted.bam"%config["named_inp"][val] for val in config["named_inp"]]),
        extname= " ".join([val for val in config["named_ext"]]),
        inpname= " ".join([val for val in config["named_inp"]])
    threads:
        10
    log:
        stdout="logs/deeptools/all_multiBAMSummary.log",
        stderr="logs/deeptools/all_multiBAMSummary.err"

    shell:
        "multiBamSummary bins -b {params.inpbams} {params.extbams} -o {output} "
        "--labels {params.inpname} {params.extname} "
        "--binSize 1000 "
        "--numberOfProcessors {threads} --extendReads "
        "--minMappingQuality 10 "
        "--samFlagInclude 66 > {log.stdout} 2> {log.stderr}"

rule deeptools_QC_corHeatmap:
    input:
        "deeptools_QC/all_coverage_matrix.npz"
    output:
        plot="deeptools_QC/all_corHeatmap.png",
        cormat="deeptools_QC/all_corvalues.txt"
    log:
        stdout="logs/deeptools/all_corHeatmap.log",
        stderr="logs/deeptools/all_corHeatmap.err"
    shell:
        "plotCorrelation --corData {input} --corMethod 'spearman' "
        "--whatToPlot 'heatmap' --plotFile {output.plot} "
        "--colorMap 'viridis' "
        "--outFileCorMatrix {output.cormat} > {log.stdout} 2> {log.stderr}"

rule deeptools_QC_plotPCA:
    input:
        "deeptools_QC/all_coverage_matrix.npz"
    output:
        plot="deeptools_QC/all_plotPCA.png",
        pcaout="deeptools_QC/all_PCA_data.txt"
    log:
        stdout="logs/deeptools/all_plotPCA.log",
        stderr="logs/deeptools/all_plotPCA.err"
    shell:
        "plotPCA --corData {input} "
        "--plotFile {output.plot} "
        "--outFileNameData {output.pcaout} > {log.stdout} 2> {log.stderr}"

rule deeptools_QC_corScatterplot:
    input:
        "deeptools_QC/all_coverage_matrix.npz"
    output:
        plot="deeptools_QC/all_corScatterplot.png"
    log:
        stdout="logs/deeptools/all_corScatterplot.log",
        stderr="logs/deeptools/all_corScatterplot.err"
    shell:
        "plotCorrelation --corData {input} --corMethod 'spearman' "
        "--whatToPlot 'scatterplot' --plotFile {output.plot} "
        "> {log.stdout} 2> {log.stderr}"

rule deeptools_QC_plotCoverage:
    input: expand("bowtie2/{sample}_sorted.bam", sample=[config["named_ext"][val] for val in config["named_ext"]]),
           expand("bowtie2/{sample}_sorted.bam.bai", sample=[config["named_ext"][val] for val in config["named_ext"]]),
           expand("bowtie2/{sample}_sorted.bam", sample=[config["named_inp"][val] for val in config["named_inp"]]),
           expand("bowtie2/{sample}_sorted.bam.bai", sample=[config["named_inp"][val] for val in config["named_inp"]]),
    output:
        plot="deeptools_QC/all_plotCoverage.png",
        counts="deeptools_QC/all_plotCoverage_count.txt"
    params:
        extbams = " ".join(["bowtie2/%s_sorted.bam"%config["named_ext"][val] for val in config["named_ext"]]),
        inpbams = " ".join(["bowtie2/%s_sorted.bam"%config["named_inp"][val] for val in config["named_inp"]]),
        extname= " ".join([val for val in config["named_ext"]]),
        inpname= " ".join([val for val in config["named_inp"]])
    threads:
        10
    log:
        stdout="logs/deeptools/all_plotCoverage.log",
        stderr="logs/deeptools/all_plotCoverage.err"

    shell:
        "plotCoverage -b {params.inpbams} {params.extbams} -o {output.plot} "
        "--labels {params.inpname} {params.extname} "
        "--numberOfProcessors {threads} --extendReads "
        "--numberOfSamples 1000 "
        "--outRawCounts {output.counts} "
        "--minMappingQuality 10 "
        "--samFlagInclude 66 > {log.stdout} 2> {log.stderr}"

rule deeptools_QC_plotCoverage_rmdups:
    input: expand("bowtie2/{sample}_sorted.bam", sample=[config["named_ext"][val] for val in config["named_ext"]]),
           expand("bowtie2/{sample}_sorted.bam.bai", sample=[config["named_ext"][val] for val in config["named_ext"]]),
           expand("bowtie2/{sample}_sorted.bam", sample=[config["named_inp"][val] for val in config["named_inp"]]),
           expand("bowtie2/{sample}_sorted.bam.bai", sample=[config["named_inp"][val] for val in config["named_inp"]]),
    output:
        plot="deeptools_QC/all_plotCoverage_rmdups.png",
        counts="deeptools_QC/all_plotCoverage_count_rmdups.txt"
    params:
        extbams = " ".join(["bowtie2/%s_sorted.bam"%config["named_ext"][val] for val in config["named_ext"]]),
        inpbams = " ".join(["bowtie2/%s_sorted.bam"%config["named_inp"][val] for val in config["named_inp"]]),
        extname= " ".join([val for val in config["named_ext"]]),
        inpname= " ".join([val for val in config["named_inp"]])
    threads:
        10
    log:
        stdout="logs/deeptools/all_plotCoverage_rmdups.log",
        stderr="logs/deeptools/all_plotCoverage_rmdups.err"
    shell:
        "plotCoverage -b {params.inpbams} {params.extbams} --plotFile {output.plot} "
        "--labels {params.inpname} {params.extname} "
        "--numberOfProcessors {threads} --extendReads --ignoreDuplicates "
        "--numberOfSamples 1000 "
        "--outRawCounts {output.counts} "
        "--minMappingQuality 10 "
        "--samFlagInclude 66 > {log.stdout} 2> {log.stderr}"

rule deeptools_QC_bamPEFragmentSize:
    input: expand("bowtie2/{sample}_sorted.bam", sample=[config["named_ext"][val] for val in config["named_ext"]]),
           expand("bowtie2/{sample}_sorted.bam.bai", sample=[config["named_ext"][val] for val in config["named_ext"]]),
           expand("bowtie2/{sample}_sorted.bam", sample=[config["named_inp"][val] for val in config["named_inp"]]),
           expand("bowtie2/{sample}_sorted.bam.bai", sample=[config["named_inp"][val] for val in config["named_inp"]]),
    output:
        plot="deeptools_QC/all_bamPEFragmentSize.png",
        fragment_hist="deeptools_QC/all_bamPEFragmentSize.txt",
        table="deeptools_QC/all_bamPEFragmentSize_table.txt"
    params:
        extbams = " ".join(["bowtie2/%s_sorted.bam"%config["named_ext"][val] for val in config["named_ext"]]),
        inpbams = " ".join(["bowtie2/%s_sorted.bam"%config["named_inp"][val] for val in config["named_inp"]]),
        extname= " ".join([val for val in config["named_ext"]]),
        inpname= " ".join([val for val in config["named_inp"]])
    threads:
        10
    log:
        stdout="logs/deeptools/all_bamPEFragmentSize.log",
        stderr="logs/deeptools/all_bamPEFragmentSize.err"
    shell:
        "bamPEFragmentSize -b {params.inpbams} {params.extbams} --histogram {output.plot} "
        "--samplesLabel {params.inpname} {params.extname} "
        "--numberOfProcessors {threads} "
        "--distanceBetweenBins 10000 "
        "--table {output.table} --outRawFragmentLengths {output.fragment_hist}  > {log.stdout} 2> {log.stderr}"

rule deeptools_QC_computeGCBias:
    input: 
        "bowtie2/{sample}_sorted.bam"
    output:
        gcfile="deeptools_QC/{sample}_GCBias_freqs.txt"
        #gcplot="deeptools_QC/{sample}_GCBias_plot.png"
    threads:
        5
    log:
        stdout="logs/deeptools/{sample}_GCBias.log",
        stderr="logs/deeptools/{sample}_GCBias.err"
    shell:
        "computeGCBias -b {input} --effectiveGenomeSize 4641652 "
        "-g /mnt/scratch/mbwolfe/genomes/U00096_3.2bit "
        "--sampleSize 50000 "
        "--numberOfProcessors {threads} "
        "--GCbiasFrequenciesFile {output.gcfile}  > {log.stdout} 2> {log.stderr}"

rule deeptools_ComputeMatrix_ref:
    input:
       lambda wildcards: ["deeptools_coverage/%s_%s.bw"%(val, wildcards.norm) for val in config["named_ext"]]
    output:
        outfile="deeptools_summary/computeMatrix_ref_{gff}_{norm}.gz",
        outmat="deeptools_summary/computeMatrix_ref_{gff}_{norm}.tab"

    wildcard_constraints:
        norm="(RPKM|CPM|BPM|RPGC|count|SES|median).*"
    threads:
        5
    log:
        stdout="logs/deeptools/computeMatrix_ref_{gff}_{norm}.log",
        stderr="logs/deeptools/computeMatrix_ref_{gff}_{norm}.err"
    params:
        bwfiles = lambda wildcards: " ".join(
                    ["deeptools_coverage/%s_%s.bw"%(val,wildcards.norm) 
                    for val in config["named_ext"]]),
        gfffile = lambda wildcards: config["gffs"][wildcards.gff],
        names = " ".join([val for val in config["named_ext"]])
    shell:
        "computeMatrix reference-point -S {params.bwfiles} -R {params.gfffile} " 
        "-a 1500 -b 1500 --referencePoint TSS --samplesLabel {params.names} "
        "--numberOfProcessors {threads} --outFileName {output.outfile} "
        "--outFileNameMatrix {output.outmat} > {log.stdout} 2> {log.stderr}"

rule deeptools_ComputeMatrix_scale:
    input:
       lambda wildcards: ["deeptools_coverage/%s_%s.bw"%(val, wildcards.norm) for val in config["named_ext"]]
    output:
        outfile="deeptools_summary/computeMatrix_scale_{gff}_{norm}.gz",
        outmat="deeptools_summary/computeMatrix_scale_{gff}_{norm}.tab"

    wildcard_constraints:
        norm="(RPKM|CPM|BPM|RPGC|count|SES|median).*"
    threads:
        5
    log:
        stdout="logs/deeptools/computeMatrix_scale_{gff}_{norm}.log",
        stderr="logs/deeptools/computeMatrix_scale_{gff}_{norm}.err"
    params:
        bwfiles = lambda wildcards: " ".join(
                    ["deeptools_coverage/%s_%s.bw"%(val,wildcards.norm) 
                    for val in config["named_ext"]]),
        gfffile = lambda wildcards: config["gffs"][wildcards.gff],
        names = " ".join([val for val in config["named_ext"]])
    shell:
        "computeMatrix scale-regions -S {params.bwfiles} -R {params.gfffile} " 
        "-a 1000 -b 1000 --samplesLabel {params.names} "
        "-m 1000 " 
        "--numberOfProcessors {threads} --outFileName {output.outfile} "
        "--outFileNameMatrix {output.outmat} > {log.stdout} 2> {log.stderr}"

rule deeptools_plotProfile_ref:
    input:
        "deeptools_summary/computeMatrix_ref_{gff}_{norm}.gz"
    output:
        outfile="deeptools_summary/plotProfile_ref_{gff}_{norm}.png",
        outmat="deeptools_summary/plotProfile_ref_{gff}_{norm}.tab"
    wildcard_constraints:
        norm="(RPKM|CPM|BPM|RPGC|count|SES|median).*"
    log:
        stdout="logs/deeptools/plotProfile_ref_{gff}_{norm}.log",
        stderr="logs/deeptools/plotProfile_ref_{gff}_{norm}.err",
    params:
        label = lambda wildcards: wildcards.gff
    shell:
        "plotProfile --matrixFile {input} --outFileName {output.outfile} "
        "--outFileNameData {output.outmat} --refPointLabel 'start' "
        "--perGroup --regionsLabel {params.label} > {log.stdout} 2> {log.stderr}"


rule deeptools_plotProfile_ref_heatmap:
    input:
        "deeptools_summary/computeMatrix_ref_{gff}_{norm}.gz"
    output:
        outfile="deeptools_summary/plotProfile_ref_{gff}_{norm}_heatmap.png",
    wildcard_constraints:
        norm="(RPKM|CPM|BPM|RPGC|count|SES|median).*"
    log:
        stdout="logs/deeptools/plotProfile_ref_{gff}_{norm}_heatmap.log",
        stderr="logs/deeptools/plotProfile_ref_{gff}_{norm}_heatmap.err",
    params:
        label = lambda wildcards: wildcards.gff
    shell:
        "plotProfile --matrixFile {input} --outFileName {output.outfile} "
        "--refPointLabel 'start' "
        "--plotHeight 10 --plotWidth 10 "
        "--perGroup --plotType heatmap --regionsLabel {params.label} > {log.stdout} 2> {log.stderr}"

rule deeptools_plotProfile_scale:
    input:
        "deeptools_summary/computeMatrix_scale_{gff}_{norm}.gz"
    output:
        outfile="deeptools_summary/plotProfile_scale_{gff}_{norm}.png",
        outmat="deeptools_summary/plotProfile_scale_{gff}_{norm}.tab"
    wildcard_constraints:
        norm="(RPKM|CPM|BPM|RPGC|count|SES|median).*"
    log:
        stdout="logs/deeptools/plotProfile_scale_{gff}_{norm}.log",
        stderr="logs/deeptools/plotProfile_scale_{gff}_{norm}.err",
    params:
        label = lambda wildcards: wildcards.gff
    shell:
        "plotProfile --matrixFile {input} --outFileName {output.outfile} "
        "--outFileNameData {output.outmat} --startLabel 'start' "
        "--endLabel 'end' "
        "--perGroup --regionsLabel {params.label} > {log.stdout} 2> {log.stderr}"

