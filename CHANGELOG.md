# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Develop

### Added
- Added a log2 strand bias calculation as a possible coverage track
- Calculation of the count of fragments for each contig
- Allow two separate cutadapt runs to occur
- Added a rule to generate annotation files easily with `run_annotations`
- Added the ability to specify relative coordinates in `bwtools_query`
- Unstranded Peak calling and motif calling
- Modeling module with DESeq2 and stranded feature quant with HTseq-count

### Changed
- Refactored `NETseq_pause_calling.py` for performance boost

### Bug fixes
- Fixed issue with `workflow/rules/coverage_and_norm.smk` extending reads for
  non-paired end samples
- Fixed issue with `multiqc` report not properly reporting pre and post trim
  fastqc reports

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

