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

- preprocessing to remove adapters and low quality sequences using [`cutadapt`](https://cutadapt.readthedocs.io/en/stable/)
  and [`trimmomatic`](http://www.usadellab.org/cms/?page=trimmomatic)
- alignment to a reference genome(s) using [`bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- read coverage calculation with a variety of normalization options within or between samples 
  with [`deeptools`](https://deeptools.readthedocs.io/en/latest/) and custom python scripts 
- peak calling using either [`macs2`](https://github.com/macs3-project/MACS) or a custom python 
  implementation of [`CMARRT`](https://github.com/mikewolfe/CMARRT_python)

In addition this pipeline performs and use [`multiqc`](https://multiqc.info/) to compile the 
following quality control metrics into an interactive html report:
- Read QC using [`fastqc`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) both before and after preprocessing
- A large number of ChIP quality control metrics calculated by [`deeptools`](https://deeptools.readthedocs.io/en/latest/)


# Installing the pipeline

This pipeline has three dependencies:
- The package manager [`miniconda`](https://docs.conda.io/en/latest/miniconda.html)
- The workflow management system [`snakemake`](https://snakemake.readthedocs.io/en/stable/index.html)
- The API for working with Portable Encaspulated Projects (PEP) [`peppy`](http://peppy.databio.org/en/latest/)

`miniconda` can be installed following the installation instructions for your system [here](https://docs.conda.io/en/latest/miniconda.html).

Once `miniconda` is installed, both `snakemake` and `peppy` can be installed in their own environment easily using
```
conda create -n ChIPseq_pipeline snakemake=5.24.2 peppy
conda activate ChIPseq_pipeline
```

Now you can pull the pipeline from github using:
```
git clone --recurse-submodules https://github.com/mikewolfe/ChIPseq_pipeline/
```
And you can change into the newly cloned `ChIPseq_pipeline` directory to test your installation with:
