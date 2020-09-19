# control of the pipeline
configfile: "config/config.yaml"
# sample metadata and information
pepfile: "pep/config.yaml"

## GLOBAL HELPER FUNCTIONS
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

def determine_extracted_samples(pep):
    samp_table = pep.sample_table
    samples = samp_table.loc[~samp_table["input_sample"].isna(), "sample_name"]
    return samples.tolist()
    

# include in several rules here
include: "workflow/rules/preprocessing.smk"
include: "workflow/rules/alignment.smk"
include: "workflow/rules/coverage_and_norm.smk"
include: "workflow/rules/quality_control.smk"



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

