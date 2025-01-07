#' ---
#' title: "Run DESeq2 for differential expression"
#' author: "Mike Wolfe"
#' output: github_document
#' ---

#' note that this is inspired from set up script
#' https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/blob/master/workflow/scripts/sleuth-init.R
#' and diffexp script
#' https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/blob/master/workflow/scripts/sleuth-diffexp.R

#' load in required libraries
suppressMessages(library(DESeq2))
suppressMessages(library("BiocParallel"))
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
metadata <- args[1]
HTseq_dir <- args[2]
outputprefix <- args[3]
full_model <- args[4]
reduced_model <- args[5]
cores <- as.integer(args[6])
spike_regions <- args[7]
if(spike_regions == "NA"){
    spike_regions <- NA
}
register(MulticoreParam(cores))
#' read in metadata table
#' make sure all strings are factors except the path and sample fields
#' TODO control ordering of factors
metadata <- read_tsv(metadata, col_names = TRUE) %>% 
                    mutate_at(vars(fileName, sample), as.character)
#' read in the full model
full_mod_form <- as.formula(full_model)
#'read the reduced model
reduced_mod_form <- as.formula(reduced_model)
# pull all possible needed variables
variables <- labels(terms(full_mod_form)) %>% strsplit('[:*]') %>% unlist()
#' make sure only samples that actually have values for everything in model are
#' included
metadata <- metadata %>% drop_na(c(sample, all_of(variables))) %>%
    select(sample, fileName, all_of(variables)) %>% as.data.frame()

#' to start we won't do any target mapping but could change that in the future
#' also this assumes that you want transcript level summaries not gene. This
#' will need to be updated in the future
glimpse(metadata)

#' create our DEseq2 matrix
dds <- DESeqDataSetFromHTSeqCount(sampleTable = metadata,
                                  directory = HTseq_dir,
                              design = full_mod_form)
#' Adjust for spike-ins if they exist
if(!is.na(spike_regions)){
    spike_regions <- read_tsv(spike_regions, col_names = c("chrm", "start", "end", "region", "value", "strand")) %>%
                            mutate(id = str_c(chrm, region, sep = "_")) %>%pull(id)
    spike_regions <- rownames(dds) %in% spike_regions
    dds <- estimateSizeFactors(dds, controlGenes = spike_regions)
}

#' fit full model
print(str_c("Fitting full model", " ", full_model))
dds <- DESeq(dds, parallel = TRUE)

#' let's get the values for each of the coefficients
covariates <- resultsNames(dds)
#' ignore the intercept
covariates <- covariates[2:length(covariates)]

#' do shrinkage based estimates of each coefficient
all_wald <- data.frame(region = rownames(dds))
for(covariate in covariates){
    print(str_c("Performing wald test for ", covariate))
    resLFC <- lfcShrink(dds, coef=covariate, type="apeglm")
    beta_col_name <- str_c("b", covariate, sep = "_")
    beta_se_col_name <- str_c(beta_col_name, "se", sep = "_")
    beta_pval_col_name <- str_c(beta_col_name, "pval", sep = "_")
    beta_qval_col_name <- str_c(beta_col_name, "qval", sep = "_")
    print(as.data.frame(resLFC) %>% as_tibble() %>% mutate(region = rownames(resLFC)) %>% select(region, everything()))
    this_wald <- as.data.frame(resLFC) %>% as_tibble() %>% mutate(region = rownames(resLFC)) %>% 
    select(region,
           !!beta_col_name := log2FoldChange,
           !!beta_se_col_name := lfcSE,
           !!beta_pval_col_name := pvalue,
           !!beta_qval_col_name := padj)
    all_wald <- inner_join(all_wald, this_wald, by = "region")
}
#' coefficients
write_tsv(all_wald, str_c(outputprefix, "_coefficients.tsv"))

#' get normalized variance stabilized counts
vsd <- vst(dds, blind = FALSE)
write_tsv(assay(vsd) %>% as_tibble() %>% 
          mutate(region = rownames(assay(vsd))) %>% 
                 select(region, everything()), str_c(outputprefix, "_scaledcount.tsv"))

#' write out estimated size factors
write_tsv(sizeFactors(dds) %>% as_tibble() %>% mutate(sample = names(sizeFactors(dds))), str_c(outputprefix, "_scalefactors.tsv"))

#' final full data structure
saveRDS(dds, str_c(outputprefix, ".rds"))


#' Do a likelihood ratio test
dds <- DESeq(dds, test = "LRT", reduced = reduced_mod_form, parallel = TRUE)
res <- results(dds)
write_tsv(res %>% as_tibble() %>%
          mutate(region = rownames(res)) %>%
          select(region, everything()), str_c(outputprefix, "_lrt.tsv"))

