# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Develop

### Added
- Pseudoalignment with Kallisto for stranded data
- Modeling with Sleuth
- Seperation of feature counting from modeling
- Stranded querys with sense and antisense capabilities `bwtools_query_stranded`
- Stranded spike-in normalization using variety of methods including Bonhoure et
  al. method and DEseq2 method
- Deduplication with Picard MarkDuplicates.
- Ability to set specific locations to zero when doing NETseq pause calling
- Give background corrected logos in NETseq

### Changes
- Update versions of some packages including ncbi-acc-download to get around
  rate filtering issue and multiqc to deal with python updates
- Enable the ability to choose which QC you want to run

### Bug fixes
- Issue with `bwtools_multiprocessing` not correctly identifying needed input
  files
- Fixed issues with sample sheets that do not have a input_sample or only have
  single end samples
- Fixed bug in NETseq pause calling that would only report last chromosome

## 0.0.3

### Added
- Added a log2 strand bias calculation as a possible coverage track
- Calculation of the count of fragments for each contig
- Allow two separate cutadapt runs to occur for sequencing kits that require
  sequential trimming
- Added a rule to generate annotation files easily with `run_annotations`
- Added the ability to specify relative coordinates in `bwtools_query`
- Unstranded Peak calling and motif calling
- Modeling module with DESeq2 and stranded feature quant with HTseq-count
- Ability to better set parameters for ChIP QC in config file

### Changed
- Refactored `NETseq_pause_calling.py` for performance boost

### Bug fixes
- Fixed issue with `workflow/rules/coverage_and_norm.smk` extending reads for
  non-paired end samples
- Fixed issue with `multiqc` report not properly reporting pre and post trim
  fastqc reports
- Fixed issue with bamPEFragmentSize would attempt to run on single end data
- Fixed issue with trimmomatic not properly running for single end data
- Fixed issue with spearman correlation calculations failing on updated coverage
  files from bwtools
- Fixed issue with NETseq pause calling failing when NaNs were present in data
  by coercing all NaNs to zeros
- Fixed issue where pulling sequences from bed files near genome boundaries
  would fail. Updated to consider all chromosomes as circular.
- Fixed an issue where NETseq logos would display information content on
  a non-standard scale. Fixed scaling to be 0 to 2.0 bits
- Fixed an issue where numpy deprecations impact deeptools quality control. Set
  a lower version of numpy in environment files for quality control

## 0.0.2

### Added
- Add ability to get combined annotation table and bed for genome annotations

### Bug fixes
- Fixed bug preventing stranded coverage under some normalization schemes

## 0.0.1
- Added ability to process coverage in a stranded manner
- Added support for NET-seq data
- Added limited support for bisulfite data


# ChIPseq pipeline change log before fork

## 0.2.8 - 2021-01-20

### Added
- Ability to collapse multiple bigwig files into one using summary stats
- Ability to perform operations between two groups of bigwig files

## 0.2.7 - 2021-01-10

### Added
- Added support for single end libraries

## 0.2.6 - 2021-12-08

### Added
- Added relative polymerase progression calculation to `bwtools_query` (http://dx.doi.org/10.1016/j.molcel.2014.02.005)
- Add the ability to smooth raw signals using convolution with a Gaussian or flat kernel. 
- Add the ability to smooth raw signals using a Savitzky-Golay filter
- Ability to calculate the region-level spearman correlations for a set of samples in `bwtools_query`
- Add calculations of the traveling ratio (10.1016/j.molcel.2008.12.021)
- Add discovery of a local max within a region

### Changed
- Output of `bwtools_query` is now compressed (gzip). File extension changed to `.tsv.gz`

## 0.2.5 - 2021-10-25

### Added
- Allow mix and match between `querysub`, `queryscale` `fixedsub`, and `fixedscale`.

## 0.2.4 - 2021-10-13

### Added
- Create a compiled html table of all mutations called across all samples in `variant_calling`

### Bug fixes
- Fixed a bug in renaming `sample/output/output.gd` files to `renamed_output/sample.gd` in `variant_calling`

## 0.2.3 - 2021-10-12

### Added
- A `variant_calling` module to call changes from the reference genome with ChIP input samples ([#2](https://github.com/mikewolfe/ChIPseq_pipeline/issues/2))
- Add option to control what fraction of a region can have NaNs in it and still report a summary value in `bwtools query`

