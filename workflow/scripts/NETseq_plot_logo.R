library(ggseqlogo)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
print(args)
d <- read_tsv(args[1], col_names = FALSE)
p <- ggseqlogo(d)

p <- p + theme_classic() 
p$scales$scales[[1]] <- scale_x_continuous(breaks = seq(1, 18, by = 1), 
                                           labels= c(seq(-13,-1,by=1), seq(1, 5, by = 1)))

ggsave(args[2], p, width = 5, height = 2)
