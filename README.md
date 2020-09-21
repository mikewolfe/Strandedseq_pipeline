# A fully portable and reproducible ChIP-seq pipeline

The goal of this project is to create a reproducible analysis pipeline for ChIP-seq data 
that can be easily deployed in any computational environment. Currently this pipeline is
aimed at processing ChIP-seq samples with the following 
characteristics:

- Paired-end short-read illumina sequencing data
- Organisms with a single chromosome

However, support for multiple chromosomes or reference genomes that include additional
contigs is a high priority development goal.

# What the pipeline does
Starting from raw fastq files this pipeline does the following:

- preprocessing to remove adapters and low quality sequences using `cutadapt` and `trimmomatic`
- alignment to a reference genome(s) using `bowtie2`
- read coverage calculation with a variety of normalization options within or between samples 
  with `deeptools` and custom python scripts 
- peak calling using either `macs2` or a custom python implementation of `CMARRT`

In addition this pipeline performs and compiles the following quality control into an interactive html report:
- Read QC using `fastqc` both before and after preprocessing
- A large number of ChIP quality control metrics calculated by `deeptools` (https://deeptools.readthedocs.io/en/develop/index.html)


# Installing the pipeline

