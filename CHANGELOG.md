# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Full support for genomes with multiple contigs/chromosomes ([#5](https://github.com/mikewolfe/ChIPseq_pipeline/issues/5))
- Full support of masking regions of the genome from analysis ([#11](https://github.com/mikewolfe/ChIPseq_pipeline/issues/11))
- Full support for normalization based on a known background region ([#6](https://github.com/mikewolfe/ChIPseq_pipeline/issues/6))
- Calculate summaries over regions or locations ([#14](https://github.com/mikewolfe/ChIPseq_pipeline/issues/14))

### Changed
- Config syntax for specifying genome inputs
- Made the pseudocount argument explicit in ratio calculations. Default before was 1.

### Bug fixes
- Fixed an issue when running `coverage_and_norm` module only ([#3](https://github.com/mikewolfe/ChIPseq_pipeline/issues/3))
- Fixed an issue with the `tbb` version being too high on some systems in the `alignment.yaml` file
- Fixed an issue where genome size was not properly read as input to Macs2 in `workflow/rules/peak_calling.smk`
