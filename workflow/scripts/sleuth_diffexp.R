#' ---
#' title: "Run sleuth for differential expression"
#' author: "Mike Wolfe"
#' output: github_document
#' ---

#' note that this is inspired from set up script
#' https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/blob/master/workflow/scripts/sleuth-init.R
#' and diffexp script
#' https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/blob/master/workflow/scripts/sleuth-diffexp.R

#' load in required libraries
suppressMessages(library(tidyverse))
suppressMessages(library(sleuth))

args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
outputprefix <- args[2]
full_model <- args[3]
reduced_model <- args[4]
cores <- args[5]
exclude <- args[6]

#' read in metadata table
#' make sure all strings are factors except the path and sample fields
#' TODO allow for samples to be excluded
#' TODO control ordering of factors
samples <- read.table(infile, stringsAsFactors = TRUE, sep = "\t",
                      header = TRUE) %>% mutate_at(vars(path, sample), as.character)
#' read in the full model
full_mod_form <- as.formula(full_model)
#'read the reduced model
reduced_mod_form <- as.formula(reduced_model)
# pull all possible needed variables
variables <- labels(terms(full_mod_form)) %>% strsplit('[:*]') %>% unlist()
#' make sure only samples that actually have values for everything in model are
#' included
samples <- samples %>% drop_na(c(sample, path, all_of(variables))) %>%
    select(sample, path, all_of(variables))

#' to start we won't do any target mapping but could change that in the future
#' also this assumes that you want transcript level summaries not gene. This
#' will need to be updated in the future
glimpse(samples)
so <- sleuth_prep(samples,
                  full_model = full_mod_form,
                  extra_bootstrap_summary = TRUE,
                  read_bootstrap_tpm = TRUE,
                  # puts everything on a log2 scale not a natural log scale
                  # the default for sleuth is log(x + 0.5)
                  transform_fun_counts = function(x) log2(x + 0.5),
                  num_cores = cores
                  )

#' fit full model
print(str_c("Fitting full model", " ", full_model))
so <- sleuth_fit(so, full_mod_form, "full")
print(str_c("Fitting reduced model", " ", reduced_model))
#' fit reduced model
so <- sleuth_fit(so, reduced_mod_form, "reduced")
#' get likelihood ratio test
print("Performing likelihood ratio test")
so <- sleuth_lrt(so, 'reduced', 'full')
all <- sleuth_results(so, "reduced:full", "lrt", show_all = TRUE,
                      pval_aggregate = FALSE) %>% arrange(pval)

#' let's get the values for each of the coefficients
covariates <- colnames(design_matrix(so, "full"))
#covariates <- covariates[covariates != "(Intercept)"]

#' following along with code from the reference above
#' to pull out all of the coefficients
all_wald <- all %>% select(target_id)
for(covariate in covariates){
    print(str_c("Performing wald test for ", covariate))
    so <- sleuth_wt(so, covariate, "full")
    beta_col_name <- str_c("b", covariate, sep = "_")
    beta_se_col_name <- str_c(beta_col_name, "se", sep = "_")
    beta_pval_col_name <- str_c(beta_col_name, "pval", sep = "_")
    beta_qval_col_name <- str_c(beta_col_name, "qval", sep = "_")

    this_wald <- sleuth_results(so, covariate, "wt", show_all = TRUE,
                               pval_aggregate = FALSE) %>%
    select(target_id = target_id,
           !!beta_col_name := b,
           !!beta_se_col_name := se_b,
           !!beta_pval_col_name := pval,
           !!beta_qval_col_name := qval)
    all_wald <- inner_join(all_wald, this_wald, by = "target_id")
}
#' likelihood ratio test
write.table(all, str_c(outputprefix, "_lrt.tsv"), col.names = TRUE, row.names = FALSE,
            quote = FALSE, sep = "\t")
#' coefficients
write.table(all_wald, str_c(outputprefix, "_coefficients.tsv"),
            col.names = TRUE, row.names = FALSE,
            quote = FALSE, sep = "\t")
#' normalized tpms
sleuth_matrix <- sleuth_to_matrix(so, 'obs_norm', 'tpm')
sleuth_matrix <- as_tibble(sleuth_matrix) %>% 
mutate(target_id = rownames(sleuth_matrix)) %>% 
select(target_id, everything())
write.table(sleuth_matrix, str_c(outputprefix, "_normed_tpm.tsv"),
            col.names = TRUE, row.names = FALSE,
            quote = FALSE, sep = "\t")

#' raw tpms
sleuth_matrix <- sleuth_to_matrix(so, 'obs_raw', 'tpm')
sleuth_matrix <- as_tibble(sleuth_matrix) %>% 
mutate(target_id = rownames(sleuth_matrix)) %>% 
select(target_id, everything())
write.table(sleuth_matrix, str_c(outputprefix, "_raw_tpm.tsv"),
            col.names = TRUE, row.names = FALSE,
            quote = FALSE, sep = "\t")

#' final full data structure
sleuth_save(so, str_c(outputprefix, ".rds"))
