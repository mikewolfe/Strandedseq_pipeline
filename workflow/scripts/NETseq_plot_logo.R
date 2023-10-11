library(ggseqlogo)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
print(args)
d <- read_tsv(args[1], col_names = c("seq"))
prior <- read_delim(args[3], col_names = c("value", "prior"), comment = "#", delim = " ")

# Manually calculate the PWM
bg_cor <- d %>% mutate(seq = str_split(seq, "", n = 18, simplify = TRUE)) %>% 
    as.matrix() %>% as_tibble() %>% pivot_longer(everything()) %>% filter(value != "N") %>%
    mutate(name = as.numeric(str_remove(name, "seq."))) %>%
    group_by(name, value) %>% summarize(n=n()) %>% 
    group_by(name) %>% mutate(freq = n/sum(n)) %>% 
    left_join(prior, by = "value") %>% 
    group_by(name) %>% 
    mutate(total_inf = sum(freq*log2(freq/prior))) %>% ungroup() %>% 
    mutate(rel_ent = pmax(freq* log2(freq/prior), 0)) %>%
    mutate(height = rel_ent) %>% select(name, value, height) %>% 
    pivot_wider(names_from = name, values_from = height) %>%
    column_to_rownames("value") %>%
    as.matrix()

p <- ggseqlogo(d, seq_type = 'dna')
p <- p + theme_classic() 
p$scales$scales[[1]] <- scale_x_continuous(breaks = seq(1, 18, by = 1), 
                                           labels= c(seq(-13,-1,by=1), seq(1, 5, by = 1)))

ggsave(str_c(args[2], ".pdf"),  p, width = 5, height = 2)

p <- ggseqlogo(bg_cor, method = "custom", seq_type = 'dna')
p <- p + theme_classic() + labs(y = "Relative Entropy")
p$scales$scales[[1]] <- scale_x_continuous(breaks = seq(1, 18, by = 1), 
                                           labels= c(seq(-13,-1,by=1), seq(1, 5, by = 1)))




ggsave(str_c(args[2], "_bg_corrected.pdf"), p, width = 5, height = 2)
