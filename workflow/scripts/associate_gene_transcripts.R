#' load in required libraries
suppressMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)

#' As input we will take a single bed file and a ptt and rnt file
#' Goal is to step through each record and determine which records in the 
#' ptt/rnt file it fully overlaps
inbed <- read_tsv(args[1], col_names = c("chrm", "start", "end", "target_id",
                                            "value", "strand"))
ptt <- read_tsv(args[2], skip = 2) %>% 
    separate(Location,into = c("Gene_start", "Gene_end"), sep = "\\.\\.") %>% 
    mutate(Gene_start = as.numeric(Gene_start), Gene_end = as.numeric(Gene_end)) %>%
    mutate(Gene_start = Gene_start - 1) %>%
    select(Gene_start, Gene_end, Strand, Gene, Synonym, Product) 
rnt <- read_tsv(args[3], skip = 2) %>%
    separate(Location,into = c("Gene_start", "Gene_end"), sep = "\\.\\.") %>% 
    mutate(Gene_start = as.numeric(Gene_start), Gene_end = as.numeric(Gene_end)) %>%
    mutate(Gene_start = Gene_start - 1) %>%
    select(Gene_start, Gene_end, Strand, Gene, Synonym, Product)
gene_info <- bind_rows(ptt, rnt)


# use a function to determine which genes overlap with which transcripts
find_overlaps <- function(t_start, t_end, t_strand, test_entries){
    out <- test_entries %>% filter(Strand == t_strand) %>%
        filter((Gene_start >= t_start) & (Gene_end <= t_end)) %>%
        select(-Strand)
    out
}

# map this over the entire set of transcripts to find matches

out <- inbed %>% 
    mutate(genes = 
               pmap(list(start, end, strand), 
                    function(x, y, z) find_overlaps(x, y, z, gene_info))) %>%
    unnest(genes, keep_empty = TRUE)

# write out the final tsv
write_tsv(out, args[4])
