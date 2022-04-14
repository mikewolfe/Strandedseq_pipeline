# NET-seq example - Larson et al. 2014 

This pipeline with example config files has everything you need to go from raw
fastqs to pause calling assuming standard illumina Tru-seq adapters and no
funkiness with the reads.  The pipeline has the flexibility to mess with funky
reads if you need it though.  

If you are running the whole pipeline, you need to download the raw fastq files 
from GEO. You can do this by running the `raw_data/fetch_data.sh` shell script
assuming that you have `fasterq-dump` installed
(https://github.com/ncbi/sra-tools/wiki/HowTo%3A-fasterq-dump).

After downloading the raw reads you can then run the NET-seq relevant portions
of the pipeline with:
`snakemake run_NETseq --use-conda --cores 10`


However, it is likely that you already have some processing pipelines
you are using so you don’t need this whole thing, just the pause calling.
The pause calling is in `workflow/scripts/NETseq_pause_calling.py`. 

You can check the help for it like so:
```
(bioinf) mwolfe6@ad.wisc.edu (Larson_NETseq) scripts $ python3 NETseq_pause_calling.py --help
usage: Python script to call pauses in NET-seq data [-h] {Genomewide,Region} ...

positional arguments:
  {Genomewide,Region}  Supported Commands
    Genomewide         Call pauses across the entire genome
    Region             Call pauses across based on defined regions

optional arguments:
  -h, --help           show this help message and exit
```

This tells you that it runs in two modes, the Genomewide mode and the Region
mode:
* Genomewide: which slides a 201 bp window over the genome and considers pauses in the context of that window
* Region-based:  which considers pauses in the context of a defined set of
                 regions. Here I used all possible ORFs in the *E coli* genome as
                 my set of regions so each pause within in ORF is considered
                 against using the whole ORF as the window. Ideally one would
                 want to pad out ORFs in a future analysis to include upstream
                 regions of the transcript involved in promoter proximal pausing

You can see the usage information for a mode like the following:

```
(bioinf) mwolfe6@ad.wisc.edu (Larson_NETseq) scripts $ python3 NETseq_pause_calling.py Genomewide --help
usage: Python script to call pauses in NET-seq data Genomewide [-h] [--p P] [--method METHOD] [--n_boot N_BOOT] [--seed SEED]
                                                               [--filters FILTERS [FILTERS ...]] [--avg_cov_cutoff AVG_COV_CUTOFF]
                                                               [--cov_cutoff COV_CUTOFF] [--zero_cutoff ZERO_CUTOFF] [--circular CIRCULAR]
                                                               [--alpha ALPHA]
                                                               infile_plus infile_minus wsize

positional arguments:
  infile_plus           bigwig file containing raw counts on the plus strand
  infile_minus          bigwig file containing raw counts on the plus strand
  wsize                 Half the window size. Adds to each side of the query bp. I.e. 100 would give you a total of a 201 bp window with the query bp
                        in the center. Default = 100

optional arguments:
  -h, --help            show this help message and exit
  --p P                 Number of processors to use. Default = 1
  --method METHOD       Must be one of Larson (gaussian), Poisson, or bootstrap. Default = Poisson
  --n_boot N_BOOT       Number of bootstraps to do for bootstrap method. Default = 10000
  --seed SEED           Random number generator starting seed for reproducibility. Default = 42
  --filters FILTERS [FILTERS ...]
                        List of filters to apply before considering windows. Possibilities include: 'zero', 'max', 'avg_cov', 'cov'. Default = None
  --avg_cov_cutoff AVG_COV_CUTOFF
                        Cutoff for average coverage filter. Default = 10
  --cov_cutoff COV_CUTOFF
                        Cutoff for center coverage filter. Default = 1
  --zero_cutoff ZERO_CUTOFF
                        Cutoff for fraction of zeros in a window. Default = .10
  --circular CIRCULAR   Consider the genome as circular? Wraps edge windows around end of genome. Default = True
  --alpha ALPHA         q-value cutoff for significance. Default = 0.05
```

Note the `--filters` option allows you to filter out windows you don’t want to
test. For example if you wanted to filter out windows that have at least 1 read
at the center and an average of 1 read across the window you would specify
`---filter avg_cov cov --cov_cutoff 1 --avg_cov_cutoff 1`. Note that if you do
not specify the `--filters` option than no filters will be applied.  

This script depends on additional custom python modules that are in that same
directory. So if you want to move the script from its location you should just
copy the whole `workflow/scripts/` directory to ensure that the `bed_utils.py`
and `bwtools.py` dependencies can be found by the script. Additionally the
script depends on the `pyBigwig` python module, the `scipy` module, and the
`numpy` module so you must be working in an environment that has those python
modules installed.  

The script takes genomewide `.bw` files as input. These can be easily generated
from your `.bam` alignment files using the `deeptools` python module. For
single-end data for RNA prepared with a dUTP type library (rf-stranded) the
commands would be:

```
bamCoverage --bam input.bam --outFileName output_plus_strand.bw --outFileFormat 'bigwig' --binSize 1 --filterRNAstrand forward --Offset 1 
bamCoverage --bam input.bam --outFileName output_minus_strand.bw --outFileFormat 'bigwig' --binSize 1 --filterRNAstrand reverse --Offset 1 
```

If you were doing paired end data you would want to additionally add a
`--samFlagInclude 67` to the command to make sure that you don’t count twice
for each read pair. Additionally if you have a fr-stranded library than the
`--filterRNAstrand` flag should just be flipped for each strand. More details
can be found in the documentation for this tool:
https://deeptools.readthedocs.io/en/stable/content/tools/bamCoverage.html 

I have implemented three different statistical methods for calling pauses:
* Larson : which considers how rare it would be to see a total count at the
           pause site given the mean and stdev of the window and the null
           hypothesis that the counts should be normally distributed if there is
           no pause.  
* Poisson: which considers how rare it would be to see a total count at the
           pause site given the average count across the window and the null
           hypothesis that the counts should be poisson distributed if there is
           no pause.  
* bootstrap : which considers how rare it would be to see a total count at the 
              pause site given a null hypothesis that the read should be evenly
              distributed across the window. This null distribution is built
              using 10,000 random draws from a multinomial distribution and the
              total count in the window. The p-value is then generated as the
              actual count + 1/ (total random draws with a value >= actual
              count) + 1

It should be noted that running the bootstrapping method genomewide with
minimal window filtering takes about 10 hours on 10 cores on a single sample,
so this thing isn’t super fast. If we need more speed I will need to implement
the algorithm in a compiled language instead.

Please let me know if you have any additional questions or concerns. Also, if
there is a feature you want (for example, filtering for windows that are only
in specific regions during the genomewide method) please let me know and I can
add it to the script.
