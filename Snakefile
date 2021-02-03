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

def determine_resolution(config):
    if "coverage" in config and "resolution" in config["coverage"]:
        resolution = config["coverage"]["resolution"]
    else:
        logger.warning(
        """
        Could not find specification for resolution in config file. I.e.

        coverage:
            resolution: N

        defaulting to 5 bp resolution
        """)
        resolution = 5
        
    return resolution

def determine_effective_genome_size_file(sample, config, pep):
    genome = lookup_sample_metadata(sample, "genome", pep)
    return "results/alignment/combine_fasta/%s/%s_mappable_size.txt"%(genome, genome)

def determine_effective_genome_size(sample, config, pep):
    infile = determine_effective_genome_size_file(sample, config, pep)
    with open(infile, mode = "r") as inf:
        size = inf.readline().rstrip()  
    return size

def determine_masked_regions_file(config, genome):
    if "masked_regions" in config["reference"][genome]:
        outfile = config["reference"][genome]["masked_regions"]
    else:
        outfile = None
    return outfile

def determine_within_normalization(config):
    if "normalization" in config and "within" in config["normalization"]:
        within = config["normalization"]["within"]
    else:
        logger.warning(
        """
        Could not find specification for normalization within a sample in config file. I.e.

        normalization:
            within: X

        defaulting to normalization by median signal in bins
        """)
        within = "median"
    return within

def determine_final_normalization(config):
    ending = "log2ratio"
    if "normalization" in config and "RobustZ" in config["normalization"]:
        RZ = config["normalization"]["RobustZ"]
        if RZ:
            ending += "RZ"
    return ending
    
RES = determine_resolution(config)
WITHIN = determine_within_normalization(config)
ENDING = determine_final_normalization(config)

# include in several rules here
include: "workflow/rules/preprocessing.smk"
include: "workflow/rules/alignment.smk"
include: "workflow/rules/coverage_and_norm.smk"
include: "workflow/rules/quality_control.smk"
include: "workflow/rules/peak_calling.smk"



## overall rules

rule run_all:
    input: 
        "results/quality_control/multiqc_report.html",
        expand("results/peak_calling/cmarrt/{sample}_{within}_{ending}.narrowPeak",\
        sample = determine_extracted_samples(pep),\
        within = WITHIN,\
        ending = ENDING),
        expand("results/peak_calling/macs2/{sample}_peaks.xls",\
        sample = determine_extracted_samples(pep))

rule clean_all:
    threads: 1
    shell:
        "rm -rf results/"    
