library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
print(args)
outfile <- args[1]
infiles <- args[2:length(args)]

read_and_combine_data <- function(infiles){
    out <- vector('list', length(infiles))
    for (num in seq_along(infiles)){
        these_data <- read_tsv(infiles[num], col_types ="ci")
        these_data <- these_data %>% mutate(sample_name = basename(infiles[num]) %>% str_replace("_frags_per_contig.tsv", ""))
        out[[num]] <- these_data
    }
    bind_rows(out)
}

d <- read_and_combine_data(infiles)
write_tsv(d, file=outfile)
